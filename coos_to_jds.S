.intel_syntax noprefix

.global coos_to_jds

/*
 * This function transforms a coos matrix to a jds matrix. It is identical to the C implementation of
 *
 * coos_to_jds_c(coos_matrix* coos, uint64_t* nnz_per_row, jds_matrix* jds)
 *                     rdi                      rsi                 rdx
 */
coos_to_jds:

	push r12
	push r13
	push r14
	push r15

/*
 * Push number of rows to stack. We will need this later. It is pushed to the stack to save registers
 * for other variables which will be accessed more often.
 */
	mov rcx, [rdi]


/*
 *		    rsp
 *			 |
 *		   	 v
 * stack := [num_rows, ...]
 */
	push rcx

/*
 * r8 := JDS permutation array address. This values in this array are initialized in utility.c with [ 0, 1, 2, ..., num_rows-1]
 */
	mov r8, [rdx + 56]

/*
 * Prepare parameters to call 'sort' while maintaining calling convention.
 * The array containing the number of non-zeros per row is sorted in descending order.
 * This is necessary to construct a correct JDS matrix.
 * This 'sort' function also sorts the permutation array to keep track of the original permutation of the rows.
 */
	push r8
	push rdx
	push rsi
	push rdi

	mov rdi, rsi
	mov rsi, rcx
	mov rdx, r8

	call sort

/*
 * Next, the triples of COOS matrix are sorted by their column index.
 * This is necessary since rows need to contain elements from left to right.
 * This step is relevant when forming the jagged diagonals of the resulting JDS matrix.
 */

	mov rdi, [rsp]

	call sort_triples_col

	pop rdi
	pop rsi
	pop rdx
	pop r8

/*		    rsp
 *			 |
 *			 v
 * stack := [num_entries, num_rows, ...]
 */
	mov rax, [rdi + 16]
	push rax

	mov [rdx + 16], rax

/*
 * Next we push the number of jagged diagonals to the stack. This number is equal to the maximum
 * number of non-zeros per row of the COOS matrix. Since the array of non-zeros per row is sorted in descending order,
 * we can simply push the first number in the array on the stack.
 *
 *		    rsp
 *			 |
 *			 v
 * stack := [number_of_jagged_diagonals, num_entries, num_rows, ...]
 */
	mov rax, [rsi]
	push rax

/*
 * Number of elements in the col_ptr array in the result JDS matrix is
 * equal to num_of_jagged_diagonals + 1
 */
	add rax, 1
	mov [rdx + 40], rax

/*
 * The first entry in the column pointer array is always 0
 */
	mov rax, [rdx + 48]
	mov qword ptr [rax], 0

/*
 * rax := values iterator
 * (used to keep track of where to write the float value in the values array in the result JDS matrix)
 */
	xor rax, rax
/*
 * rcx := jagged diagonals iterator (used to loop through the jagged diagonals).
 * Its value is also needed to when filling up the column pointer array.
 * rcx is .Lloop1's iterator
 */
	xor rcx, rcx
.Lloop1:
/*
 * r9 := nonzeros_per_row array iterator.
 * (used to decrement the number of nonzeros per row as we fill up 1 jagged diagonal.
 * It's also used to access the permutation array of the result JDS matrix)
 * r9 is .Lloop2's iterator
 */
	xor r9, r9
.Lloop2:
/*
 * r10 := skip
 * (used when iterating through the triples as we fill up 1 jagged diagonal. We need to skip over a certain number triples
 * for a fixed row before choosing a triple to be added to the result)
 */
	xor r10, r10
/*
 * r11 := triples iterator (used to iterate through all the COOS triples of the parameter matrix as we construct the result JDS matrix)
 * r11 is .Lloop3's iterator
 */
	xor r11, r11
.Lloop3:
	mov r12, [rdi + 24]

	push rax
	push rdx

	mov rax, 24

	mul r11
	mov r13, rax

	pop rdx
	pop rax

	add r12, r13
	mov r13, [r12]

/*
 * Check if the row index of the current triple equals the row in the permutation array indexed by r9.
 */
	cmp r13, [r8 + r9 * 8]
	jne .Lelse2
/*
 * Check if we have skipped over enough elements. For the n-th jagged diagonal,
 * we need to skip over n-1 triples of a certain row before using the values of that triple to fill our JDS matrix.
 */
	cmp r10, rcx
	jne .Lelse

/*
 * If we've indeed skipped over enough triples, we need to enter its value and column index in the respective arrays in the JDS result.
 */
	mov r14d, [r12 + 16]
	mov r13, [rdx + 24]
	mov [r13 + rax * 4], r14d

	mov r14, [r12 + 8]
	mov r13, [rdx + 32]
	mov [r13 + rax * 8], r14
/*
 * Since we just added elements to the values and column index array of the result JDS matrix respectively,
 * we need to increment the values iterator (rax). (Note that both values and col_idx array of the JDS format have the same size).
 */
	add rax, 1
	jmp .Lbreak

.Lelse:
/*
 * Since we found a triple with a row index equal to the current row in the permutation array, we increment the skip variable (r10). This will help us check
 * how many triples with a certain row index need to be skipped over before we fill up the JDS result matrix
 */
	add r10, 1
.Lelse2:
/*
 * Increment .Lloop3 iterator and check loop condition
 */
	add r11, 1
	cmp r11, [rsp + 8]
	jb .Lloop3
.Lbreak:
	sub qword ptr [rsi + r9 * 8], 1

/*
 * Increment .Lloop2 iterator and check loop condition
 */
	add r9, 1
	cmp r9, [rsp + 16]
	jae .Lcontinue
	cmp qword ptr [rsi + r9 * 8], 0
	jz .Lcontinue
	jmp .Lloop2
.Lcontinue:
	mov r12, [rdx + 48]
	mov [r12 + rcx * 8 + 8], rax

/*
 * Increment .Lloop1 iterator and check loop condition
 */
	add rcx, 1
	cmp rcx, [rsp]
	jb .Lloop1
.Lret:
/*
 * Return stack pointer to its original position
 */
	add rsp, 24
	pop r15
	pop r14
	pop r13
	pop r12
	ret
