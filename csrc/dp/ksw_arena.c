/* Bump arena for ksw2 (see ksw_arena.h). Each block has a 16-byte header holding its size,
 * for realloc; payloads are 16-byte aligned. */
#include "ksw_arena.h"

#include <stdlib.h>
#include <string.h>

#define KSW_ARENA_ALIGN ((size_t)16)
#define KSW_ARENA_HEADER ((size_t)16)
/* Ceiling on one arena chunk: minimap2's 256 MB kalloc block size. A full-band gap fill's
 * traceback matrix, (q+t) x min(q, 2w+1) bytes, fits below the max_sw_mat cap; a larger
 * request goes to malloc instead of pinning its peak in every worker. */
#define KSW_ARENA_CAPACITY_LIMIT ((size_t)256 << 20)

struct ksw_arena_s {
  unsigned char* block; /* the malloc'd chunk (may be NULL) */
  unsigned char* base;  /* 16-byte aligned start inside `block` */
  size_t capacity;      /* usable bytes from `base` */
  size_t offset;        /* bump cursor */
  size_t wanted;        /* bytes requested since the last reset */
};

static size_t ksw_arena_round(size_t size) {
  return (size + (KSW_ARENA_ALIGN - 1)) & ~(KSW_ARENA_ALIGN - 1);
}

ksw_arena_t* ksw_arena_create(void) {
  return (ksw_arena_t*)calloc(1, sizeof(ksw_arena_t));
}

void ksw_arena_destroy(ksw_arena_t* arena) {
  if (!arena)
    return;
  free(arena->block);
  free(arena);
}

size_t ksw_arena_capacity(const ksw_arena_t* arena) {
  return arena ? arena->capacity : 0;
}

void ksw_arena_reset(ksw_arena_t* arena) {
  if (!arena)
    return;
  if (arena->wanted > arena->capacity &&
      arena->capacity < KSW_ARENA_CAPACITY_LIMIT) {
    size_t want = arena->wanted > KSW_ARENA_CAPACITY_LIMIT
                      ? KSW_ARENA_CAPACITY_LIMIT
                      : arena->wanted;
    unsigned char* grown = (unsigned char*)malloc(want + KSW_ARENA_ALIGN);
    if (grown) { /* grow-only; a failed growth simply keeps the old chunk */
      free(arena->block);
      arena->block = grown;
      arena->base = (unsigned char*)(((size_t)grown + (KSW_ARENA_ALIGN - 1)) &
                                     ~(KSW_ARENA_ALIGN - 1));
      arena->capacity = want;
    }
  }
  arena->offset = 0;
  arena->wanted = 0;
}

static int ksw_arena_owns(const ksw_arena_t* arena, const void* ptr) {
  const unsigned char* p = (const unsigned char*)ptr;
  return arena->base != 0 && p >= arena->base + KSW_ARENA_HEADER &&
         p < arena->base + arena->capacity;
}

void* ksw_arena_malloc(void* km, size_t size) {
  ksw_arena_t* arena = (ksw_arena_t*)km;
  size_t need;
  if (!arena)
    return malloc(size);
  need = KSW_ARENA_HEADER + ksw_arena_round(size);
  arena->wanted += need;
  if (arena->offset + need <= arena->capacity) {
    unsigned char* block = arena->base + arena->offset;
    memcpy(block, &size, sizeof(size));
    arena->offset += need;
    return block + KSW_ARENA_HEADER;
  }
  return malloc(size); /* does not fit: system allocator */
}

void* ksw_arena_calloc(void* km, size_t count, size_t size) {
  size_t total;
  void* block;
  if (!km)
    return calloc(count, size);
  if (size != 0 && count > (size_t)-1 / size)
    return 0;
  total = count * size;
  block = ksw_arena_malloc(km, total);
  if (block)
    memset(block, 0, total);
  return block;
}

void* ksw_arena_realloc(void* km, void* ptr, size_t size) {
  ksw_arena_t* arena = (ksw_arena_t*)km;
  if (arena && ptr && ksw_arena_owns(arena, ptr)) {
    size_t previous;
    void* block;
    memcpy(&previous, (unsigned char*)ptr - KSW_ARENA_HEADER, sizeof(previous));
    block = ksw_arena_malloc(km, size);
    if (block && previous)
      memcpy(block, ptr, previous < size ? previous : size);
    return block;
  }
  return realloc(ptr, size);
}

void ksw_arena_free(void* km, void* ptr) {
  ksw_arena_t* arena = (ksw_arena_t*)km;
  if (!ptr)
    return;
  if (arena && ksw_arena_owns(arena, ptr))
    return; /* ksw_arena_reset owns it */
  free(ptr);
}
