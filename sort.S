.intel_syntax noprefix

.global sort

/*
  This function receives 2 arrays as parameters. The first array is the array to be sorted and
  the second array simply contains the numbers 0,1,2, ..., (size of 1st array - 1).

  The function sorts the first array in descending order and swaps elements in the second
  array everytime it swaps elements in the first array.
  This allows us to keep track of the permutation of the elements in the first array as they get
  shuffled due to the sorting. The permutation is then found in the second array.

  The sorting algorithm used is insertion sort.
								^^^^^^^^^

  This function is identical to the C implementation of:
          sort_c(uint64_t* array_1, uint64_t size, uint64_t* array_2)
*/
sort:
//An array of size 1 (or 0) doesn't need to be sorted
    cmp rsi, 1
    jbe .Lret

//rax := outer loop iterator
    mov rax, 1
.Lloop:
//rcx := inner loop iterator
    xor rcx, rcx
.LloopInner:

    mov r8, [rdi + rcx * 8]
//Compare the 2 ints in the first array
    cmp r8, [rdi + rax * 8]
    jae .Lcont

//Swap elements in the first array
    mov r9, [rdi + rax * 8]
    mov [rdi + rax * 8], r8
    mov [rdi + rcx * 8], r9
//Swap elements in the second array
    mov r8, [rdx + rcx * 8]
    mov r9, [rdx + rax * 8]
    mov [rdx + rax * 8], r8
    mov [rdx + rcx * 8], r9
.Lcont:
//Inner loop condition and iterator incrementation
    inc rcx
    cmp rcx, rax
    jb .LloopInner

//Outer loop condition and iterator incrementation
    inc rax
    cmp rax, rsi
    jb .Lloop
.Lret:
    ret
