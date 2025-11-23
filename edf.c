#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define MAX_TASKS 64
#define MAX_NAME 64

typedef struct {
    char name[MAX_NAME];
    double period;
    double wcet[4];
} Task;

typedef struct {
    int task_idx;
    long job_no;
    double release;
    double deadline;
    double remaining;
} Job;

typedef struct Node {
    Job job;
    struct Node *next;
} Node;

typedef struct {
    int n;
    double horizon;
    double P_active_mW[4];
    double P_idle_mW;
    Task tasks[MAX_TASKS];
} System;

static const int FREQ_MHZ[4] = {1188, 918, 648, 384};

static int dcmp(double a, double b) {
    const double eps = 1e-9;
    if (a < b - eps) return -1;
    if (a > b + eps) return 1;
    return 0;
}

static void rq_insert(Node **head, Job j) {
    Node *n = (Node*)malloc(sizeof(Node));
    n->job = j; n->next = NULL;
    if (!*head) { *head = n; return; }
    Node *p = NULL, *c = *head;
    while (c) {
        int cmp = dcmp(j.deadline, c->job.deadline);
        if (cmp < 0 || (cmp == 0 && j.task_idx < c->job.task_idx)) break;
        p = c; c = c->next;
    }
    if (!p) { n->next = *head; *head = n; }
    else { n->next = p->next; p->next = n; }
}

static bool rq_pop_front(Node **head, Job *out) {
    if (!*head) return false;
    Node *n = *head; *out = n->job; *head = n->next; free(n); return true;
}

static void rq_free(Node **head) {
    Node *c = *head; while (c) { Node *n = c->next; free(c); c = n; } *head = NULL;
}

static void simulate_edf(const System *sys) {
    Node *rq = NULL;
    double t = 0.0;
    double horizon = sys->horizon;
    long k[MAX_TASKS] = {0};
    double next_rel[MAX_TASKS] = {0};

    for (int i = 0; i < sys->n; ++i) {
        Job j = {i, 0, 0.0, sys->tasks[i].period, sys->tasks[i].wcet[0]};
        rq_insert(&rq, j);
        k[i] = 1;
        next_rel[i] = sys->tasks[i].period;
    }

    double totalE = 0, idle = 0;
    bool have_cur = false;
    Job cur = {0};

    while (t < horizon) {
        double next_release = horizon + 1e9;
        for (int i = 0; i < sys->n; ++i)
            if (next_rel[i] < next_release) next_release = next_rel[i];

        if (!have_cur && !rq) {
            double until = fmin(next_release, horizon);
            double dt = until - t;
            if (dt > 0) {
                double e = (sys->P_idle_mW/1000.0)*dt;
                printf("%.0f IDLE IDLE %.0f %.3fJ\n", t, dt, e);
                totalE += e; idle += dt; t += dt;
            }
            continue;
        }

        for (int i = 0; i < sys->n; ++i) {
            while (dcmp(next_rel[i], t) == 0) {
                Job j = {i, k[i], t, t + sys->tasks[i].period, sys->tasks[i].wcet[0]};
                rq_insert(&rq, j);
                k[i]++; next_rel[i] += sys->tasks[i].period;
            }
        }

        if (!have_cur) { rq_pop_front(&rq, &cur); have_cur = true; }

        double next_event = fmin(t + cur.remaining, next_release);
        double dt = fmin(next_event, horizon) - t;
        double e = (sys->P_active_mW[0]/1000.0)*dt;
        printf("%.0f %s 1188 %.0f %.3fJ\n", t, sys->tasks[cur.task_idx].name, dt, e);
        totalE += e;
        cur.remaining -= dt; t += dt;

        if (t > cur.deadline && cur.remaining > 1e-9) {
            fprintf(stderr, "Deadline miss %s\n", sys->tasks[cur.task_idx].name);
            rq_free(&rq); return;
        }
        if (cur.remaining <= 1e-9) have_cur = false;
    }

    printf("TOTAL_ENERGY %.3fJ\nPCT_IDLE %.2f%%\nEXEC_TIME %.0fs\n",
           totalE, (idle*100.0/sys->horizon), sys->horizon);
}

int main(int argc, char **argv) {
    if (argc < 2) { printf("Usage: %s <input.txt>\n", argv[0]); return 1; }
    FILE *fp = fopen(argv[1], "r");
    if (!fp) { perror("fopen"); return 1; }
    System sys = {0};
    fscanf(fp, "%d %lf %lf %lf %lf %lf %lf", &sys.n, &sys.horizon,
           &sys.P_active_mW[0], &sys.P_active_mW[1], &sys.P_active_mW[2],
           &sys.P_active_mW[3], &sys.P_idle_mW);
    for (int i = 0; i < sys.n; ++i)
        fscanf(fp, "%s %lf %lf %lf %lf %lf",
               sys.tasks[i].name, &sys.tasks[i].period,
               &sys.tasks[i].wcet[0], &sys.tasks[i].wcet[1],
               &sys.tasks[i].wcet[2], &sys.tasks[i].wcet[3]);
    fclose(fp);
    simulate_edf(&sys);
}
