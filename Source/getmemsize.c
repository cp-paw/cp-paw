/* -- AIX/6000 System monitor 
**
**     fps.c
**
** Copyright (c) 1991-1995 Jussi Maki, All Rights Reserved.
** Copyright (c) 1993-1996 Marcel Mol, All Rights Reserved.
** NON-COMMERCIAL USE ALLOWED. YOU ARE FREE TO DISTRIBUTE
** THIS PROGRAM AND MODIFY IT AS LONG AS YOU KEEP ORIGINAL
** COPYRIGHTS.
*/

/* Fast ps like utility "fps" -- Jussi Maki, Center for Scientific Computing,
 * Finland, jmaki@csc.fi,  February 21, 1992 - July 15, 1996
 * compile this program in AIX with eg. command:
 *   cc -o fps fps.c
 */

#include <stdio.h>
#include <procinfo.h>
#define NPROCS 10000

/* system call to get process table */
extern getproc(struct procinfo *procinfo, int nproc, int sizproc);
    /* procinfo:   pointer to array of procinfo struct */
    /* nproc:      number of user procinfo struct */
    /* sizproc:    size of expected procinfo structure */

/* system call to get user-area variables according to process */
extern getuser(struct procinfo *procinfo, int plen, void *user, int ulen);
    /* procinfo:   ptr to array of procinfo struct */
    /* plen:       size of expected procinfo struct */
    /* user:       ptr to array of userinfo struct, OR user */
    /* ulen:       size of expected userinfo struct */

/* memory allocation doesn't really happen until it is touched
 * so these large arrays don't waste memory in AIX 3
 */
struct procinfo proc1[NPROCS];
struct userinfo user1[NPROCS];

/*
 * Function Delcarations
 */
double getmemsize(int pid);
double mem_of_proc(pid_t pid, struct procinfo *proc,
                                 struct userinfo *user, int nproc);
int    getprocessinfo(struct procinfo *procinfo, struct userinfo *userinfo);

double getmemsize(int pid)
{
  int nproc1; 
  nproc1 = getprocessinfo(proc1, user1);
  return mem_of_proc(pid, proc1, user1, nproc1);
}

double mem_of_proc(pid_t pid, struct procinfo *proc, struct userinfo *user,
               int nproc)
{
    int i;
    double value=0.0E0;

    for (i = 0; i < nproc; i++)
        if (proc[i].pi_pid == pid)
            break;
    if (i >= nproc) {
        fprintf(stderr, "No such process.\n"
                        "Maybe NPROCS is to low, then recompile me...\n");
        exit(0);
    }
    value = user[i].ui_tsize / 1024.0 + user[i].ui_dvm * 4;
    return value;
} /* mem_of_proc */

int getprocessinfo(struct procinfo *procinfo, struct userinfo *userinfo)
{
    int nproc;
    char *swappername = "swapper";
    int i;

    /* get the whole process table */
    nproc = getproc(procinfo, NPROCS, sizeof(struct procinfo));
    for (i = 0; i < nproc; i++) /* get each user-area entrys by process */
        getuser(&procinfo[i], sizeof(struct procinfo),
	        &userinfo[i], sizeof(struct userinfo));
    strcpy(userinfo[0].ui_comm, swappername); /* 1st process is always pid 0 */

    return nproc;

} /* getprocessinfo */
