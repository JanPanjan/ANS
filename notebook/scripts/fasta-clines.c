#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MAX_LINE_LEN 100

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: ./fasta-clines.out <fasta-file>");
        exit(1);
    }

    FILE *file = fopen(argv[1], "r");

    if (file == NULL) {
        perror("Error opening file.");
        exit(1);
    }

    int n = 0;               // line counter
    char line[MAX_LINE_LEN]; // line buffer

    // for measuring time more precisely than time_t
    struct timeval start_time;
    struct timeval end_time;
    double final_time = 0;

    printf("Processing file...\n");
    gettimeofday(&start_time, NULL);

    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] == '>') {
            n++;
        }
    }

    gettimeofday(&end_time, NULL);
    // first term is difference in seconds
    // second term is diference in microseconds converted to seconds
    final_time = (double)(end_time.tv_sec - start_time.tv_sec) +
        (double)(end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("sequences in file: %d\n", n);
    printf("time: %f sec\n", final_time);

    return 0;
}
