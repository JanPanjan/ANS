#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LEN 1000

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: ./gtf-genes.out <gtf-file>");
        exit(1);
    }

    FILE *file = fopen(argv[1], "r");

    if (file == NULL) {
        perror("Error opening file.");
        exit(1);
    }

    int g = 0;               // gene counter
    char line[MAX_LINE_LEN]; // line buffer

    printf("Processing file...\n");

    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] == '#') { // skip comments at the beginning of the file
            continue;
        }
        line[strcspn(line, "\n")] = 0;

        int col_num = 0; // column counter
        char *token;     // this will point to a string that represents an attribute

        token = strtok(line, "\t"); // take a token out of the string

        while (token != NULL) {
            col_num++;

            if (col_num == 3) {
                // returns 0 if it matches, 1 if not
                g += !strcmp(token, "gene");
                break;
            }

            token = strtok(NULL, "\t");   // take the next token
        }
    }

    fclose(file);

    printf("Genes in file: %d\n", g);

    return 0;
}