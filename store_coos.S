.intel_syntax noprefix

.global store_coos
/*
 * This function loops through the triples of the parameter matrix and adds the parameter value to the
 * correct row and column. If there is no such triple with matching row and column indices, we simply append
 * a new COOS triple at the end of the result array (we also increment the number of entries in the array).
 *
 * void store_coos(uint64_t row, uint64_t col, float value, void* result)
 *                     rdi             rsi            xmm0        rdx
 */

store_coos:
/*
 * Size of 1 COOS triple in bytes (including padding at the end of each element)
 */
    mov r9, 24
/*
 * Move number of entries in rcx. This will be needed for our loop condition.
 */
    mov rcx, [rdx + 16]
/*
 * rax := Loop iterator
 */
    xor rax, rax
.Lloop:
/*
 * r10 := The address of array of triples in the parameter COOS matrix
 */
    mov r10, [rdx + 24]
    push rax
    push rdx
/*
 * Get array offset in bytes and add it to the array address (r10)
 */
    mul r9
    add r10, rax
    pop rdx
    pop rax

/*
 * Compare loop iterator (rax) with number of triples (rcx). If they're equal that means we have iterated through
 * all of the triples of the parameter matrix and did not find a triple with matching row and column indices. Hence, we
 * have to add a new triple at the end of the number of entries
 */
    cmp rcx, rax
    je .Lnew_entry

/*
 * Compare row index of current triple with parameter row (rdi). If they're not equal, move on to the next triple.
 */
    cmp [r10], rdi
    jne .Lcontinue

/*
 * Compare column index of current triple with parameter column (rsi). If they're not equal, move on to the next triple.
 */
    cmp [r10 + 8], rsi
    jne .Lcontinue

/*
 * Since the we found a triple with matching row and column indices, we add the paramter value (xmm0)
 * to the COOS triple's value and return.
 */
    addss xmm0, [r10 + 16]
    movss [r10 + 16], xmm0
    ret
.Lcontinue:
    add rax, 1
    jmp .Lloop

.Lnew_entry:
/*
 * Increase number of COOS entries since we didn't find an existing entry with the same row and column indices
 */
    add qword ptr [rdx + 16], 1
/*
 * Add a new COOS entry with (row, col, val) at the end of the parameter matrix's triples array
 */
    mov [r10], rdi
    mov [r10 + 8], rsi
    movss [r10 + 16], xmm0
    ret
