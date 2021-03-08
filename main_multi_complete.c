#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <mpi.h>
#include <omp.h>

#define MAX_KEY_SIZE  8
#define MAX_WORD_SIZE 20
#define MAX_HEX_LENGTH 16

struct stat dic_sb;
struct stat enc_sb;
typedef struct pair
{
    char key[MAX_HEX_LENGTH + 1];
    int num_of_matches;
} Pair;
typedef struct hash_node
{
    char ROOT;
    char* word;
    int word_length;
    struct hash_node* next_word;
    struct hash_node* next_letter;
} Hash_node;

int processKey(char *key);
int processThreadKey(char *key, char *threadedKeyBytes);
void encode_file(int nBytes, FILE* output_file);
void encode_file_to_mem_threaded(int nBytes, char* encode_dest, char* threadedKeyBytes);
int generate_Next_Key(char* input_curr_key);
long POW(int n, int exponent);

int check_key_matches_HASH(char* key, Hash_node* hash_entry, char* encode_dest, char* threadedKeyBytes);
Hash_node* map_dict_from_mem_to_hash(char* dictionary_pointer, int dictionary_length);
int find_word_in_hash_dictionary(char* word, int word_len, Hash_node* dict_start);
void insert_to_hash(Hash_node* node, char* new_word, int new_word_len);
void free_hash_table(Hash_node* table_start);

// argv[1] is the number of BITS in the encryption KEY
// argv[2] is the path to the encrypted file
// argv[3] optional, file containing expected unencrypted words.

char keyBytes[MAX_KEY_SIZE];
char *current_Key = NULL;
char* best_key = NULL;
int best_key_hits = 0;
int num_of_hex_letters = 0;
// Default dictionary option
char* dict_path = "/usr/share/dict/american-english";
// File descriptors
int enc_fd = -1;
int dic_fd = -1;
// File lengths
int enc_file_len = 0;
int dic_file_len = 0;
// Pointers to mapped files.
char* enc_file_in_memory = NULL; // Encrypted file
char* dict_in_memory = NULL; // Dictionary file
char* temp_file_in_memory = NULL; // Temporary file to compare with dictionary
// File structs



void main(int argc, char* argv[])
{   
    omp_set_num_threads(omp_get_max_threads());

    int key_bit_length = atoi(argv[1]);
    num_of_hex_letters = key_bit_length / 4;

    // Initialize MPI
    MPI_Datatype PairMPI;
    MPI_Datatype types[2] = {MPI_CHAR, MPI_INT};
    int block_lengths[2] = {MAX_HEX_LENGTH + 1, 1}; // for creating datatype
    MPI_Aint disp[2];
    Pair myPair;
    int my_rank,num_procs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

    disp[0] = offsetof(Pair, key);
    disp[1] = offsetof(Pair, num_of_matches);

    MPI_Type_create_struct(2, block_lengths, disp, types, &PairMPI);
    MPI_Type_commit(&PairMPI);

    Pair* pair_arr = (Pair*)calloc(num_procs, sizeof(Pair));
    if (pair_arr == NULL)
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    else
    {
        int i;
        for(i = 0; i < num_procs; i++)
            pair_arr[i].key[num_of_hex_letters] = '\0';
    }

    if (my_rank == 0) // Master
    {   // Initiallize files
        printf("num of proccesses: %d\n", num_procs);
        printf("Num of arguments: %d\nKey length: %s\n", argc, argv[1]);

        enc_fd = open(argv[2], O_RDONLY, S_IRUSR | S_IWUSR);
        if (argc > 3)
            dic_fd = open(argv[3], O_RDONLY, S_IRUSR | S_IWUSR);
        
        else
            dic_fd = open(dict_path, O_RDONLY, S_IRUSR | S_IWUSR);

        if (fstat(dic_fd, &dic_sb) == -1)
            perror("coulnd't get dictionary+ file size.\n");
        if(fstat(enc_fd, &enc_sb) == -1)
            perror("coulnd't get encrypted file size.\n");

        enc_file_len = enc_sb.st_size;
        dic_file_len = dic_sb.st_size;
        printf("Encrypted file length %d\n", enc_file_len);
        enc_file_in_memory = mmap(NULL, enc_sb.st_size, PROT_READ, MAP_PRIVATE, enc_fd, 0);
        dict_in_memory = mmap(NULL, dic_sb.st_size, PROT_READ, MAP_PRIVATE, dic_fd, 0);
    }

    // Make sure all processes have the files lengths
    MPI_Bcast(&enc_file_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dic_file_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Wait until all of them received the broadcast
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (my_rank != 0) // Master already has them mapped
    {
        enc_file_in_memory = malloc(sizeof(char) * enc_file_len);
        dict_in_memory = malloc(sizeof(char) * dic_file_len);
    }
    MPI_Barrier(MPI_COMM_WORLD); // wait until all processes allocated memory
    // Send everyone the mapped files
    MPI_Bcast(enc_file_in_memory, enc_file_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(dict_in_memory, dic_file_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes
    Hash_node* hash_table_dictionary = map_dict_from_mem_to_hash(dict_in_memory, dic_file_len);


    current_Key = malloc(sizeof(char) * (num_of_hex_letters + 1));
    
    
    int i;
    // create clean array.
    for (i = 0; i < num_of_hex_letters; i++)
    {
        current_Key[i] = '0';
    }
    current_Key[num_of_hex_letters] = '\0';
    int nbytes = 0;
    char c;
    // Base 16 = 16 options (0,F), exponent = number of available "slots" for chars.
    long num_of_possible_keys = POW(16, num_of_hex_letters);
    char* best_key = malloc(sizeof(char) * (num_of_hex_letters + 1));
    best_key[num_of_hex_letters] = '\0';
    best_key_hits = 0;
    if (my_rank == 0)    
        printf("Number of possible keys: %ld\n", num_of_possible_keys);  
    myPair.key[num_of_hex_letters] = '\0';
    MPI_Barrier(MPI_COMM_WORLD);    

    clock_t begin;
    if (my_rank == 0)
        begin = clock();
    ///////////////
    // Heavy work//
    ///////////////
    
    if (dict_in_memory != NULL && enc_file_in_memory != NULL)
    {
        int max_thread_num = omp_get_max_threads();
        int step_size = max_thread_num * num_procs; // Number of keys to "skip"
        Pair* threads_best_pair_array = (Pair*)calloc(sizeof(Pair), max_thread_num);

        #pragma omp parallel shared(threads_best_pair_array, max_thread_num, hash_table_dictionary, enc_file_in_memory)
        {   
            int my_tid = omp_get_thread_num();
            int iterator, p, completion_flag = 0;
            char* threadedKeyBytes = malloc(sizeof(char) * MAX_KEY_SIZE);
            char* decoded_file_dest = malloc(sizeof(char) * enc_file_len);
            char* best_threaded_key = malloc(sizeof(char) * (num_of_hex_letters + 1));
            char* threaded_key = malloc(sizeof(char) * (num_of_hex_letters + 1));
            int best_threaded_key_matches, threaded_key_matches;
            threaded_key[num_of_hex_letters] = '\0';
            best_threaded_key[num_of_hex_letters] = '\0';
            threads_best_pair_array[my_tid].key[num_of_hex_letters] = '\0';

            //i = my_rank + num_procs * my_tid; // Starting positions i = key decimal value.
            // Every step is num_procs * max_num_threads
            
            strcpy(threaded_key, current_Key); // Everyone has their own key.

            for (iterator = 0; iterator < my_rank + num_procs * my_tid; iterator++)
            {
                if(!generate_Next_Key(threaded_key))
                {
                    printf("triggered completion flag");
                    completion_flag = 1;
                    break;
                }
            }

            if (iterator < num_of_possible_keys && completion_flag == 0)
            {
                // Initiate default best pair to initial key.
                strcpy(best_threaded_key, threaded_key);
                best_threaded_key_matches = check_key_matches_HASH(best_threaded_key, hash_table_dictionary, decoded_file_dest, threadedKeyBytes);
            }

            while(iterator < num_of_possible_keys && completion_flag == 0)
            {
                // Advance key to next possible key in series.
                for (p = 0; p < step_size && iterator < num_of_possible_keys; p++, iterator++)
                {
                    if(!generate_Next_Key(threaded_key))
                    {
                        completion_flag = 1;
                        break;
                    }
                }
                if (iterator >= num_of_possible_keys || completion_flag == 1)
                {
                    break;
                }
                
                threaded_key_matches = check_key_matches_HASH(threaded_key, hash_table_dictionary, decoded_file_dest, threadedKeyBytes);

                if (threaded_key_matches > best_threaded_key_matches)
                {
                    best_threaded_key_matches = threaded_key_matches;
                    strcpy(best_threaded_key, threaded_key);
                }
            }

            threads_best_pair_array[my_tid].num_of_matches = best_threaded_key_matches;
            strcpy(threads_best_pair_array[my_tid].key, best_threaded_key);


            if(decoded_file_dest)
                free(decoded_file_dest);
            if (threaded_key)
                free(threaded_key);
            if (best_threaded_key)
                free(best_threaded_key);
            if (threadedKeyBytes)
                free(threadedKeyBytes);

            #pragma omp barrier
        }
        best_key_hits = threads_best_pair_array[0].num_of_matches;
        strcpy(best_key, threads_best_pair_array[0].key);
        for (i = 1; i < max_thread_num; i++)
        {
            if (threads_best_pair_array[i].num_of_matches > best_key_hits)
            {
                best_key_hits = threads_best_pair_array[i].num_of_matches;
                strcpy(best_key, threads_best_pair_array[i].key);
            }
        }

        strcpy(myPair.key, best_key);
        myPair.num_of_matches = best_key_hits;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&myPair, 1, PairMPI, pair_arr, 1, PairMPI, 0, MPI_COMM_WORLD);
    
        free(threads_best_pair_array);

        if (my_rank == 0)
        {
            clock_t end = clock();
            // Find best key from gathered keys
            for ( i = 0; i < num_procs; i++)
            {
                if (best_key_hits < pair_arr[i].num_of_matches)
                    {
                        best_key_hits = pair_arr[i].num_of_matches;
                        strncpy(best_key, pair_arr[i].key, num_of_hex_letters);
                    }
            }
            if (best_key_hits > 0)
            {
                printf("Key: %s, Matches: %d\n", best_key, best_key_hits);
                printf("Execution time: %f\n", (double)(end-begin) / CLOCKS_PER_SEC);
                // allocate memory for the temp file which will hold the "decoded" file to print
                temp_file_in_memory = malloc(sizeof(char) * enc_file_len);
                nbytes = processKey(best_key);
                encode_file(nbytes, NULL);
                fwrite(temp_file_in_memory, sizeof(char), enc_file_len, stdout);
            }   
        }
    }

    

    


    // Free allocations.
    if (temp_file_in_memory != NULL)
        free(temp_file_in_memory);
    if (hash_table_dictionary != NULL)
        free_hash_table(hash_table_dictionary);

    if (my_rank == 0)
    {
        munmap(enc_file_in_memory, enc_file_len);
        munmap(dict_in_memory, dic_file_len);
    }
    else
    {
        free(enc_file_in_memory);
        free(dict_in_memory);
    }
    
    free(current_Key);
    free(best_key);
    MPI_Finalize();

}



// convert  hexadecimal digit to number between 0 and 15
// examples:  hexint('3') returns 3,  hexint('A') returns 10, hex2int('B') returns 11
// hexa letters may be lowercase or uppercase (e.g. 'a' or 'A')
int hex2int(char h) {
    h = toupper(h); // if h is a lowercase letter then convert it to uppercase

    if (h >= '0' && h <= '9')
        return h - '0';
    else if (h >= 'A' && h <= 'F')
        return h - 'A' + 10;
    else { 
        fprintf(stderr, "key should contain hexa digits\n");
        exit(1);
    }
    return 0;
}


/**
  This function stores the bytes of the key in the array 'keyBytes'. This will simplify the 
  access to the bytes later on.

   example:  Assume the key has 4 bytes:   <byte 0> <byte 1> <byte 2> <byte 3>

   Then the contents of 'KeyBytes'  will be:
   char keyBytes[] = { <byte 0>, <byte 1>, <byte 2>, <byte 3> };  
   
   byte0  will be applied  to the first character of the plaintext,
   <byte 1> will be applied to the second character of the plaintext and so on.

   The number of bytes in the key is returned.
   (the number of hexa digits in the key must be even so that the key has an integer number of bytes)

**/

int processThreadKey(char *key, char *threadedKeyBytes)
{
    int n = strlen(key);
    if (n%2 || n/2 > MAX_KEY_SIZE) {
        fprintf(stderr, "key must have even number of bytes. Number of bytes \
                should not exceed %d\n", MAX_KEY_SIZE);
        exit(1);
    }

    for(int i = 0; i < n; i += 2) {
        threadedKeyBytes[i/2] = (hex2int(key[i]) << 4) | hex2int(key[i+1]);
    }
    return n/2;
}

int processKey(char *key) {
    int n = strlen(key);
    if (n%2 || n/2 > MAX_KEY_SIZE) {
        fprintf(stderr, "key must have even number of bytes. Number of bytes \
                should not exceed %d\n", MAX_KEY_SIZE);
        exit(1);
    }

    for(int i = 0; i < n; i += 2) {
        keyBytes[i/2] = (hex2int(key[i]) << 4) | hex2int(key[i+1]);
    }
    return n/2;
}

void encode_file_to_mem_threaded(int nBytes, char* encode_dest, char* threadedKeyBytes)
{
    int i, j;
    if (encode_dest != NULL)
    {
        for (i = 0, j = 0; j < enc_file_len ; j++) {
            encode_dest[j] = enc_file_in_memory[j] ^ threadedKeyBytes[i];
            i++;
            if (i >= nBytes)
                i = 0;
        }
    }
}
void encode_file(int nBytes, FILE* output_file)
{
    int i, j;
    if (output_file == NULL)
    {
        for (i = 0, j = 0; j < enc_file_len ; j++) {
            temp_file_in_memory[j] = enc_file_in_memory[j] ^ keyBytes[i];
            i++;
            if (i >= nBytes)
                i = 0;
        }

    }
    else
    {
        for (i = 0, j = 0; j < enc_file_len ; j++) {
        fputc(enc_file_in_memory[j] ^ keyBytes[i], output_file); // ^  is the xor operator)
        i++;
        if (i >= nBytes)
            i = 0;
        }
    }
}

int generate_Next_Key(char* input_curr_key)
{
    if (input_curr_key == NULL)
        return 0;
    else
    {
        int i;
        for (i = 1; i <= num_of_hex_letters; i++)
        {
            if (input_curr_key[num_of_hex_letters-i] >= 'F')
            {
                input_curr_key[num_of_hex_letters-i] = '0';
                continue;
            }
            else if (input_curr_key[num_of_hex_letters-i] >= '9' &&  input_curr_key[num_of_hex_letters-i] < 'A')
            {
                input_curr_key[num_of_hex_letters-i] = 'A';
                break;
            }
            else
            {
                input_curr_key[num_of_hex_letters-i] += 1;
                break;
            }
        }
        int check_reset = 1;
        for (i = 0; i < num_of_hex_letters; i++)
        {
            if (input_curr_key[i] != '0')
                check_reset = 0;
        }
        if (check_reset == 1)
            return 0;
        else
            return 1;
        
    }
    
}

long POW(int n, int exponent)
{
    long result = n;
    int i;
    for (i = 0; i < exponent-1; i++)
        result *= n;
    return result;
}

int check_key_matches_HASH(char* key, Hash_node* hash_entry, char* encode_dest, char* threadedKeyBytes)
{
    int nbytes = processThreadKey(key, threadedKeyBytes);
    encode_file_to_mem_threaded(nbytes, encode_dest, threadedKeyBytes);
    int j = 0, k = 0;
    char str_buffer[MAX_WORD_SIZE], c;
    int word_complete_flag = 0;
    int temp_word_count = 0;
    while (k < enc_file_len)
    {
        c = encode_dest[k];
        if (j < MAX_WORD_SIZE && word_complete_flag == 0)
        {
            if (c != ' ' && c != '\n' && c != '\0')
            {
                str_buffer[j] = c;
                j++;
            }
            else
            {
                str_buffer[j] = '\0';
                if (strlen(str_buffer) > 0)
                    word_complete_flag = 1;
                else // invalid word.
                    j = 0;
                
            }
        }
        if (word_complete_flag)
        {
            if(find_word_in_hash_dictionary(str_buffer, j, hash_entry)) // HEAVY FUNCTION
            {
                 temp_word_count += 1;
                //successful_byte_read += strlen(str_buffer);
            }
            // start next word.
            j = 0;
            word_complete_flag = 0;
        }
        k++;
    }

    return temp_word_count;
}
Hash_node* map_dict_from_mem_to_hash(char* dictionary_pointer, int dictionary_length)
{
    // Initialize dictionary hash table
    Hash_node* table_start = malloc(sizeof(Hash_node));
    if (table_start == NULL)
        return NULL; // Failed to allocate memory
    Hash_node* current_node = table_start;
    int i;
    for (i = 0; i < 26; i++)
    {
        if (current_node != NULL)
        {
            current_node->ROOT = 65 + i; // Set inital ROOTS to A, B, C,  etc.
            current_node->word_length = 0; // Word length = 0 as it's the ROOT.
            current_node->word = NULL;
            current_node->next_word = NULL;
            if(i < 25)
                current_node->next_letter = malloc(sizeof(Hash_node));
            else
            {
                current_node->next_letter = NULL;
            }
            current_node = current_node->next_letter;
        }
    }

    // Mapping dictionary from memory to hash table
    current_node = table_start; // Start from the beggining.
    int curr_dic_pos = 0;
    int dic_word_len = 0;
    int k, insert_flag = 0;
    char dic_char;
    while (curr_dic_pos < dictionary_length)
    {
        insert_flag = 0;

        // find length of word in dictionary
        for (k = 0; (dic_char = dictionary_pointer[curr_dic_pos + k]) != '\n' && dic_char != '\0'; k++) 
        {
            dic_word_len++;
        }

        int letter_index = 0;
        while(letter_index < 52 && insert_flag == 0 && curr_dic_pos < dictionary_length)  // Go over all letters to see if it fits any.
        {
            if (toupper(dict_in_memory[curr_dic_pos]) == current_node->ROOT)
            {   
                insert_flag = 1;
                if (dic_word_len > 0)
                {
                    insert_to_hash(current_node, dict_in_memory+curr_dic_pos, dic_word_len);
                }
            }
            else
            {   
                if (current_node->next_letter == NULL && letter_index < 52)
                {
                    current_node = table_start;
                }
                else if (current_node->next_letter == NULL) // Reached end of Alphabet
                {
                    curr_dic_pos = dictionary_length;
                    break;
                }
                current_node = current_node->next_letter;
                letter_index++;
            }
        }
        curr_dic_pos += (dic_word_len +1);
        dic_word_len = 0;
    }
    return table_start;

}
int find_word_in_hash_dictionary(char* word, int word_len, Hash_node* dict_start)
{
    if (word_len == 0 || word == NULL || dict_start == NULL) // invalid parameters
        return 0;
    
    int flag = 1, i, j, k, letter_index = ((dict_start->ROOT) - 65);
    char dic_char, current_letter;
    Hash_node* current_node = dict_start;
    
    if (toupper(word[0]) >= 65 && toupper(word[0]) <= 90) // check if first letter actually in ABC.
    {
        while(current_node->ROOT != toupper(word[0]) && letter_index < 26) // Skip to the hash table entry of relevant letter.
        {
            letter_index++;
            current_node = current_node->next_letter;
        }

        if (!current_node) // Couldn't find an entry root with same first letter as input word -> Not in dictionary.
        {
            return 0;
        }   
            
        else
        {
            while(current_node != NULL)
            {
                flag = 0;
                if (word_len == (current_node->word_length))
                {
                    flag = 1;
                    for (k = 0; k < word_len; k++)
                    {
                        if(toupper(word[k]) != toupper(current_node->word[k])) // if a single character is different
                        {
                            flag = 0;
                            break;
                        }
                    }
                    
                }
                if (flag == 1) // Flag managed to remain complete => Word found
                    return 1;
                
                current_node = current_node->next_word; // Next word in hash table of same root.
            }
            // Failed to find word.
            
        }
        
    }
    return 0;
}
void insert_to_hash(Hash_node* node, char* new_word, int new_word_len)
{
    if (node != NULL)
    {
        while(node->next_word != NULL) // Reach the last word of current ROOT
            node = node->next_word;
        node->next_word = malloc(sizeof(Hash_node)); // Allocate memory for new word hash entry
        node->next_word->word = malloc(sizeof(char) * (new_word_len + 1)); // Allocate memory for new word
        node->next_word->word[new_word_len] = '\0'; // Assign null terminator at end.
        node->next_word->ROOT = node->ROOT;
        node->next_word->word_length = new_word_len;
        node->next_word->next_letter = node->next_letter;
        strncpy(node->next_word->word, new_word, new_word_len);    
        node->next_word->next_word = NULL;
    }
}
void free_hash_table(Hash_node* table_start)
{
    Hash_node* current_node = table_start, *temp;
    Hash_node* next_lett = table_start->next_letter;
    while(current_node != NULL)
    {
        if (current_node->word != NULL)
            free(current_node->word);
        if (current_node->next_word != NULL)
        {
            temp = current_node;
            current_node = current_node->next_word;
            free(temp);
        }
        else
        {
            temp = current_node;
            current_node = current_node->next_letter;
            if (current_node != NULL)
            {
                next_lett = temp->next_letter;
            }
            free(temp);
        }
        
    }
}