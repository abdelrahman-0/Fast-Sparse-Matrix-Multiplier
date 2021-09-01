.intel_syntax noprefix

.global matr_mult_jds

matr_mult_jds:

/*
 * void matr_mult_jds_c(const void* a, const void* b, void* result)
 *							rdi			rsi				rdx
 *
 *
 * This functions compares all nonzero elements of the second matrix against all nonzero elements of the first matrix and looks for matching triples. The function stores the intermediate result as a COOS matrix.
 * This will then be transformed to a JDS matrix.
 */

/*
 * push callee-saved registers to stack
 */
	push rbx
	push rbp
	push r12
	push r13
	push r14
	push r15

	push rdx
	push rsi
	push rdi

	mov rdi, 32

	call malloc

/*
 * check if returned pointer is not NULL
 */
	cmp rax, 0
	jz .Lhandle_null_ptr
/*
 * initialize number of entries in intermediate COOS matrix to 0.
 */
	mov qword ptr [rax + 16], 0
/*
 * set number of rows
 */
	mov rcx, [rsp]
	mov rcx, [rcx]
	mov [rax], rcx
/*
 * set number of columns
 */
	mov rcx, [rsp + 8]
	mov rcx, [rcx + 8]
	mov [rax + 8], rcx
/*
 * push COOS matrix pointer to the stack
 */
	push rax
/*
 * calculate size of total COOS triples in bytes and allocate enough memory for the array of triples on the COOS matrix
 */
	mov rcx, [rsp + 8]
	mov rax, [rcx + 16]

	mov rcx, [rsp + 16]
	mov rcx, [rcx + 16]
	mul rcx
	mov rdi, rax

	push rdx

	mov r8, [rsp + 8]
	mov rax, [r8]
	mov r8, [r8 + 8]
	mul r8

	pop rdx

	cmp rdi, rax
	jbe .Lmin
	mov rdi, rax
.Lmin:
	mov rsi, 24
	call calloc
	cmp rax, 0
	jz .Lhandle_null_ptr
	mov rcx, [rsp]
	mov [rcx + 24], rax

	mov rdi, [rsp]
	mov rdi, [rdi]
	mov rsi, 8

	call calloc
	cmp rax, 0
	jz .Lhandle_null_ptr
	mov rbx, rax

	pop rax
	pop rdi
	pop rsi
	pop rdx

/*
 * r15 := iterator that iterates through each Jagged Diagonal (JD) of second matrix (.Lloop1 iterator).
 */
	xor r15, r15
.Lloop1:
/*
 * rcx := position of the current element in the current JD of the second matrix
 *
 *		 (from 0 to (length of current JD - 1)).
 */
	xor rcx, rcx
/*
 * r8 := iterates through each element of each JD of the second matrix (.Lloop2 iterator).
 */

	mov r8, [rsi + 48]
	mov r8, [r8 + r15 * 8]
.Lloop2:
/*
 * r14 :=  column index of current non zero element of second matrix. Its value corresponds to the column index
 * of where the products will be accumulated (summed up) in the result matrix.
 */
	mov r14, [rsi + 32]
	mov r14, [r14 + r8 * 8]
/*
 * r9 := iterator that iterates through each JD of the first matrix (.Lloop3 iterator).
 */
	xor r9, r9
.Lloop3:
/*
 * r10 := position of the current element in the current jagged diagonal of the first matrix.
 *		 (from 0 to (length of current JD - 1)).
 */
	xor r10, r10
/*
 * r11 := iterates through each element of each JD of the first matrix. (.Lloop4 iterator)
 */
	mov r11, [rdi + 48]
	mov r11, [r11 + r9 * 8]
.Lloop4:
/*
 * r12 := column index of current non-zero element of first matrix.
 */
	mov r12, [rdi + 32]
	mov r12, [r12 + r11 * 8]
/*
 * r13 := permutation (i.e. original row index) of current non-zero element of second matrix.
 */
	mov r13, [rsi + 56]
	mov r13, [r13 + rcx * 8]
/*
 * check and see wether the two non-zero elements of both matrices should be multiplied with each other.
 */
	cmp r12, r13
	jne .Lelse
/*
 * multiply current non-zero elements of both matrices together.
 */
	mov r12, [rdi + 24]
	movss xmm0, [r12 + r11 * 4]
	mov r13, [rsi + 24]
	movss xmm1, [r13 + r8 * 4]
	mulss xmm0, xmm1

	push rdi
	push rsi
	push rdx
	push r8
	push r9
	push r10
	push r11
	push rcx
	push rax

/*
 * prepare parameters and call store_coos
 */
	mov r13, [rdi + 56]
	mov r13, [r13 + r10 * 8]

	mov rdi, r13
	mov r12, [rax + 16]
	mov rsi, r14
	mov rdx, rax

	call store_coos

	pop rax
	pop rcx
	pop r11
	pop r10
	pop r9
	pop r8
	pop rdx
	pop rsi
	pop rdi

	cmp r12, [rax + 16]
	je .Lelse
	add qword ptr [rbx + r13 * 8], 1

.Lelse:
/*
 * increment position of current non-zero element of first matrix in its corresponding jagged diagonal
 */
	add r10, 1
/*
 * .Lloop4 condition (loop up to current JD's length (of first matrix))
 */
	add r11, 1
	mov r12, [rdi + 48]
	mov r12, [r12 + r9 * 8 + 8]
	cmp r11, r12
	jb .Lloop4
/*
 * .Lloop3 condition (loop up to number of available JDs of first matrix)
 */
	add r9, 1
	mov r12, [rdi + 40]
	sub r12, 1
	cmp r9, r12
	jb .Lloop3
/*
 * increment position of current non-zero element of second matrix in its corresponding JD
 */
	add rcx, 1
/*
 *.Lloop2 condition (loop up to current JD's length (of second matrix))
 */
	add r8, 1
	mov r12, [rsi + 48]
	mov r12, [r12 + r15 * 8 + 8]
	cmp r8, r12
	jb .Lloop2
/*
 *.Lloop1 condition (loop up to number of available JDs of second matrix)
 */
	add r15, 1
	mov r12, [rsi + 40]
	sub r12, 1
	cmp r15, r12
	jb .Lloop1
.Lret:
/*
 * pop callee-saved registers
 */

	push rdx
	push rax

	mov rdi, rax
	call eliminate_zeros

	mov rax, [rsp]
	mov rdx, [rsp + 8]

	cmp qword ptr [rax], 0
	jz .Lhandle_null_ptr
	mov rdi, rax
	mov rsi, rbx
	call coos_to_jds

    mov rdi, [rsp]
    mov rdi, [rdi + 24]
    call free

    mov rdi, [rsp]
    call free

    mov rdi, rbx
    call free

    add rsp, 16
	pop r15
	pop r14
	pop	r13
	pop r12
	pop rbp
	pop rbx

	ret
/*
 * If malloc/calloc returns a NULL ptr, we set the number of rows in the resulting jds matrix to 0 and return,
 */
.Lhandle_null_ptr:
	mov rdx, [rsp + 16]
	mov qword ptr [rdx], 0
	add rsp, 72
	ret
