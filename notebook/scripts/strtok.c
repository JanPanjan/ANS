#include <stdio.h>
#include <string.h> // Required for strtok() and strcpy() (or strdup() if preferred)

int main() {
    // The string we want to tokenize.
    // Let's use a sentence with spaces and commas as potential delimiters.
    char original_string[] = "Hello, world! This is a test string.";

    // strtok modifies the string it works on.
    // If you need to keep the original string intact, make a copy.
    char string_to_tokenize[100]; // Make sure this buffer is large enough
    strcpy(string_to_tokenize, original_string);

    // Define the delimiter(s).
    // We can specify multiple delimiters. For example, " ,.!" would tokenize by space, comma, period, or exclamation mark.
    // For this first example, let's just use space as a delimiter.
    const char *delimiters = " "; // Delimiter is a space character

    char *token; // This pointer will point to the beginning of each token

    printf("Original string: \"%s\"\n", original_string);
    printf("String to tokenize: \"%s\"\n", string_to_tokenize);
    printf("Tokenizing with delimiter: \"%s\"\n\n", delimiters);

    // 1. First call to strtok():
    //    - The first argument is the string to be tokenized.
    //    - The second argument is the set of delimiter characters.
    token = strtok(string_to_tokenize, delimiters);

    int token_count = 1;
    // Loop while strtok() keeps finding tokens
    while (token != NULL) {
        printf("Token %d: %s\n", token_count, token);

        // 2. Subsequent calls to strtok():
        //    - Pass NULL as the first argument. strtok() remembers its previous position
        //      in the string it was tokenizing from the last call.
        //    - The second argument is still the set of delimiter characters.
        token = strtok(NULL, delimiters);
        token_count++;
    }

    printf("\nAfter tokenization, the 'string_to_tokenize' looks like this internally (with nulls):\n");
    for (int i = 0; i < strlen(original_string) + 1 && i < sizeof(string_to_tokenize) ; i++) {
        if (string_to_tokenize[i] == '\0') {
            printf("\\0"); // Print \0 for null terminators
        } else {
            printf("%c", string_to_tokenize[i]);
        }
    }
    printf("\n");
    printf("Notice how the spaces were replaced by null terminators in 'string_to_tokenize'.\n");
    printf("The 'original_string' remains unchanged: \"%s\"\n", original_string);


    // --- Example 2: Using multiple delimiters (e.g., space and comma) ---
    printf("\n--- Example 2: Multiple Delimiters ---\n");
    char another_string[] = "apple,banana orange;grape";
    char string_to_tokenize_2[100];
    strcpy(string_to_tokenize_2, another_string);

    const char *multi_delimiters = " ,;"; // Delimiters are space, comma, and semicolon

    printf("Original string: \"%s\"\n", another_string);
    printf("Tokenizing with delimiters: \"%s\"\n\n", multi_delimiters);

    token = strtok(string_to_tokenize_2, multi_delimiters);
    token_count = 1;
    while (token != NULL) {
        printf("Token %d: %s\n", token_count, token);
        token = strtok(NULL, multi_delimiters);
        token_count++;
    }
    printf("After tokenization 2, 'string_to_tokenize_2' is: %s (first part, others separated by \\0)\n", string_to_tokenize_2);


    // --- Example 3: Tokenizing your tab-delimited file line ---
    printf("\n--- Example 3: Tab-Delimited String ---\n");
    char file_line[] = "column1\tcolumn2\tgene_data\tcolumn4"; // Simulating a line from your file
    char line_copy[100];
    strcpy(line_copy, file_line);

    const char *tab_delimiter = "\t";

    printf("Original line: \"%s\"\n", file_line);
    printf("Tokenizing with delimiter: \"\\t\" (tab)\n\n");

    char *col_token;
    int col_num = 0;

    col_token = strtok(line_copy, tab_delimiter);
    while(col_token != NULL) {
        col_num++;
        printf("Column %d: %s\n", col_num, col_token);
        if (col_num == 3) {
            printf("   -> This is the third column: %s\n", col_token);
            // Here you would do: if (strcmp(col_token, "your_target_word") == 0) { ... }
        }
        col_token = strtok(NULL, tab_delimiter);
    }

    printf("--------------------------------------------\n");
    printf("--------------------------------------------\n");
    printf("--------------------------------------------\n");

    printf("%d\n", strcmp("gene", "gene"));
    printf("%d\n", strcmp("transcript", "gene"));

    return 0;
}