.intel_syntax noprefix

.global matr_mult_coos

/*
 * matr_mult_coos(const void* a, const void* b, void* result)
 *						rdi			rsi				rdx
 */
matr_mult_coos:
/*
 * Set number of rows and columns on the resulting matrix.
 * number of rows of resulting matrix = number of rows of first matrix
 */
	mov rax, [rdi]
	mov [rdx], rax
/*
 * number of columns of resulting matrix = number of columns of second matrix
 */
	mov rax, [rsi + 8]
	mov [rdx + 8], rax
/*
 * number of COOS triples of result matrix is initialized with 0 (this is only relevant for the storage function)
 */
	mov qword ptr [rdx + 16], 0
/*
 * iterator for outer loop
 */
	xor rax, rax
/*
 * loop conditions for both loops
 */
	push [rdi + 16]
	push [rsi + 16]
/*
 * loop through each COOS triple of first matrix
 */
.Lloop1:
/*
 * save column index
 */
	mov rcx, [rdi + 24]

/*
 * multiply the outer loop's iterator (rax) by the size of one coos_triple element (including padding) i.e by 24
 */
	mov r11, rax
	shl r11, 3
	shl rax, 4
	add rcx, r11
	add rcx, rax
	shr rax, 4
	mov r8, [rcx + 8]

/*
 * iterator for inner loop
 */
	xor r9, r9
/*
 * for each coos triple in the first matrix, loop through each COOS triple of second matrix
 */
.Lloop2:
/*
 * multiply the inner loop's iterator (r9) by the size of one coos_triple element (including padding) i.e 24
 */
	mov r10, [rsi + 24]

	mov r11, r9
	shl r11, 3
	shl r9, 4
	add r10, r11
	add r10, r9
	shr r9, 4

/*
 * compare column index (stored in rcx) with row index of current COOS triple of second matrix
 */
	cmp r8, [r10]
	jne .Lcontinue
/*
 * multiply both entries
 */
	movss xmm0, [rcx + 16]
	movss xmm1, [r10 + 16]
	mulss xmm0, xmm1

/*
 * push caller-saved registers (which are used) to stack
 */
	push rdi
	push rsi
	push rdx
	push rax
	push rcx
	push r8
	push r9
	push r10
/*
 * call function to store product at the correct row and column.
 * Paramters of store_coos are:
 *	- rdi = row index
 *	- rsi = column index
 *	- xmm0 = product (already in xmm0)
 *	- rdx = result pointer (already in rdx)
 */
	mov rdi, [rcx]
	mov rsi, [r10 + 8]
	call store_coos
/*
 * pop caller-saved registers from stack
 */
	pop r10
	pop r9
	pop r8
	pop rcx
	pop rax
	pop rdx
	pop rsi
	pop rdi
/*
 * [rcx] row idx of outer, [r10 + 8] col idx of inner
 */
.Lcontinue:
/*
 * increment inner loop's iterator
 */
	add r9, 1
	cmp r9, [rsp]
	jb .Lloop2
/*
 * increment outer loop's iterator
 */
	add rax, 1
/*
 * check outer loop's condition
 */
	cmp rax, [rsp + 8]
	jb .Lloop1
	add rsp, 16

	mov rdi, rdx

	call eliminate_zeros

	ret
