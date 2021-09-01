.intel_syntax noprefix

.global eliminate_zeros

/*
 * void eliminate_zeros(coos_matrix)

 * This function eliminates zeros triple in the parameter matrix. These zero triples could arise from equal-valued
 * positive and negative partial products being added to the same position in the matrix which then cancel out. The entry might
 * be zero, which means it should not be included in the resulting COOS matrix.

 * The main idea here is to :
 * 	1. Create a temporary array of COOS triples of size equal to the number of entries (triples) in the parameter matrix.
 * 	2. Iterate through the COOS triples of the parameter matrix and only add the triples with non-zero values to the temporary array
 * 	3. If a triple with a value of zero is found, the number of entries is then decremented.
 * 	4. The parameter matrix array pointer is freed and reassigned to the temporary array..
 */

eliminate_zeros:
    push rdi
/*
 * rax := size of 1 COOS triple. The size of 2 uint64_t and 1 float value including memory padding is 24
 */
    mov rax, 24
    mul qword ptr [rdi + 16]
/*
 * Parameter to malloc is the size to be allocated
 */
    mov rdi, rax
    call malloc

    pop rdi

/*
 * rdx := length of array of triples (number of entries). We will need this for our loop condition
 */
	mov rdx, [rdi + 16]

/*
 *Check if pointer returned from malloc is not NULL
 */
    cmp rax, 0
    jz .Lhandle_null_ptr

/*
 * Load constant 0.0 in xmm1. The values of the triples will compared with this register
 */
    xorps xmm1, xmm1

/*
 * Pointer to coos_entries of the paramater matrix
 */
    mov r9, [rdi + 24]

/*
 * The following reigsters (r11 and rsi) are byte iterators of the coos_entries array and the temporary array.
 * They are incremented by 24 after each loop.
 */

/*
 * r11 := byte iterator of original matrix
 */
    xor r11, r11
/*
 * rsi := byte iterator in temp coos_entries
 */
    xor rsi, rsi

/*
 * Iterate through all triples of parameter matrix
 */

/*
 *rcx := triples iterator
 */
	xor rcx, rcx

/*
 * If num entries of parameter matrix is 0, then go to .Lret
 */
    cmp rcx, rdx
    je .Lret
.Lloop:
/*
 * Compare value of triple to 0.0
 */
    ucomiss xmm1, dword ptr [r9 + r11 + 16]
/*
 * Zero valued triples are not added to the temporary array
 */
    jz .L_is_a_zero
/*
 * Increase number of entries by 1. (This will be decreased by 1 shortly, so the net effect is 0)
 */
    add qword ptr [rdi + 16], 1

/*
 * Move row idx to current triple in temporary array
 */
    mov r10, qword ptr [r9 + r11]
    mov qword ptr [rax + rsi], r10

/*
 * Move col idx to current triple in temporary array
 */
    mov r10, qword ptr [r9 + r11 + 8]
    mov qword ptr [rax + rsi + 8], r10

/*
 * Move value to current triple in temporary array
 */
    movss xmm0, [r9 + r11 + 16]
    movss [rax + rsi + 16], xmm0

/*
 * Increment temporary triples array's byte iterator by 24
 */
    add rsi, 24
.L_is_a_zero:
/*
 * Decrement the number of entries by 1
 */
    sub qword ptr [rdi + 16], 1
/*
 * Increment original triples array's byte iterator by 24
 */
    add r11, 24

/*
 * Increment iterator of triples
 */
    add rcx, 1

/*
 * Loop condition (iterator <= num_entries_original ?)
 */
    cmp rcx, rdx
    jb .Lloop

/*
 * Free original triples array and replace it with the temporary array
 */
    push rax
    push rdi
    mov rdi, [rdi + 24]
    call free
    pop rdi
    pop rax
    mov [rdi + 24], rax

.Lret:
    ret

/*
 * In case malloc fails, the number of rows of the paramter matrix is set to 0.
 * This is then checked outside the function and dealt with appropriately
 */
.Lhandle_null_ptr:
    mov qword ptr [rdi], 0
    ret

