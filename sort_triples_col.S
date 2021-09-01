.intel_syntax noprefix
.global sort_triples_col
.text
/*
 * This function sorts the triples of a COOS matrix by their column index. This is needed in the COOS
 * to JDS transformation. Sorting the triples by column allows us to build the JDS format correctly
 * since the elements of 1 row are then arranged from left to right.
 *
 * The sorting algorithm used is insertion sort.

 * This function is identical to the C implementation sort_triples_col_c(coos_matrix* coos).
 */

sort_triples_col:
/*
 * Calling convention
 */
    push r12

    cmp qword ptr [rdi + 16], 1
    jbe .Lret

    mov rsi, [rdi + 24]

/*
 * Two byte iterators
 */
    mov r8, 24

    mov rax, 1
.Lloop:
    xor r9, r9
    xor rcx, rcx
.LloopInner:
    mov r10, [rsi + r8 + 8]
    mov r11, [rsi + r9 + 8]
/*
 * compare arr[j] with arr[i]
 */
    cmp r11, r10
    jbe .Lcont
/*
 * swap arr[i], arr[j]
 */

/*
 * swap rows of both triples
 */
    mov rdx, [rsi + r8]
    mov r12, [rsi + r9]
    mov [rsi + r8], r12
    mov [rsi + r9], rdx
/*
 * swap cols of both triples
 */
    mov rdx, [rsi + r8 + 8]
    mov r12, [rsi + r9 + 8]
    mov [rsi + r8 + 8], r12
    mov [rsi + r9 + 8], rdx
/*
 * swap values of both triples
 */
    mov rdx, [rsi + r8 + 16]
    mov r12, [rsi + r9 + 16]
    mov [rsi + r8 + 16], r12
    mov [rsi + r9 + 16], rdx
.Lcont:
    add r9, 24
    inc rcx
    cmp rcx, rax
    jb .LloopInner

    add r8, 24
    inc rax
    cmp rax, [rdi + 16]
    jb .Lloop
.Lret:
    pop r12
    ret
