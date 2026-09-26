/* kalloc-compatible bump arena for ksw2's km argument, one per worker thread and reset once
 * per kernel call. It grows (never shrinks) to the previous call's high-water mark, up to a
 * fixed ceiling; a request that does not fit goes to the system allocator. km == NULL uses
 * the system allocator throughout. A pointer the arena does not own always goes to the
 * system allocator, which keeps ez->cigar (grown by krealloc from NULL, and outliving the
 * kernel call) on the heap. */
#ifndef KSW_ARENA_H_
#define KSW_ARENA_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ksw_arena_s ksw_arena_t;

/* Create/destroy one arena. A NULL return is not an error: it degrades to the
 * system allocator through the km == NULL path. */
ksw_arena_t* ksw_arena_create(void);
void ksw_arena_destroy(ksw_arena_t* arena);

/* Release every arena block at once and size the arena for the next call. */
void ksw_arena_reset(ksw_arena_t* arena);

/* Bytes in the arena's current chunk. */
size_t ksw_arena_capacity(const ksw_arena_t* arena);

/* kalloc-compatible entry points behind ksw2.h's kmalloc/kcalloc/... macros. */
void* ksw_arena_malloc(void* km, size_t size);
void* ksw_arena_calloc(void* km, size_t count, size_t size);
void* ksw_arena_realloc(void* km, void* ptr, size_t size);
void ksw_arena_free(void* km, void* ptr);

#ifdef __cplusplus
}
#endif

#endif
