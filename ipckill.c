#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>

int main(int argc, char* argv[])
{
  key_t key;
  int semid;
  /*  union semun arg;*/
  
  /*  if ((key = ftok("semdemo.c", 'J')) == -1) {
    perror("ftok");
    exit(1);
    }*/
  
  /* grab the semaphore set created by seminit.c: */
  /*if ((semid = semget(key, 1, 0)) == -1) {
    perror("semget");
    exit(1);
    }*/
  
  /* remove it: */
  /*  if (semctl(6160511, 0, IPC_RMID,arg) == -1) {*/

  /*  printf("argument is %s", argv[1]);*/

  int id;
  int tt;

  sscanf(argv[1], "%d", &tt);
  sscanf(argv[2], "%d", &id);
 
  if (tt == 0) {
    /* remove semaphore */
    if (semctl(id, 0, IPC_RMID) == -1) {
      perror("semctl");
      exit(1);
    }
  } else {
    /* remove shmem */
    if (shmctl(id, 0, IPC_RMID) == -1) {
      /*  if (semctl(id, 0, IPC_RMID) == -1) {*/
      perror("semctl");
      exit(1);
    }
  }
  return 0;
}

