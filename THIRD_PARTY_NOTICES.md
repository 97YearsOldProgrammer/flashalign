# Third-party notices

## minimap2 behavioral provenance

FlashAlign's native ordinary long-read RNA splice controller is an
independently structured C++ implementation of a deliberately narrow behavior
contract. Its filtering, packet-boundary, Z-drop replay, split, and CIGAR
compatibility behavior was derived in part by studying minimap2 2.30-r1287
(`align.c`) and checking against an externally built oracle. No minimap2 mapping
controller source, private header, data structure, or runtime component is
included in the production implementation.

minimap2 is distributed under the MIT License:

> Copyright (c) 2018- Dana-Farber Cancer Institute
> Copyright (c) 2017-2018 Broad Institute, Inc.
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in
> all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

## KSW2 SIMD kernels

The low-level KSW2 SIMD kernels under `csrc/dp/` are separately retained and
credited MIT-licensed implementations by Heng Li and contributors. Each source
file carries its upstream copyright and license. They are used through
FlashAlign's narrow DP interfaces and are not a mapping controller.

## klib kseq.h

`csrc/io/vendor/kseq.h` is kseq.h from klib by Attractive Chaos
(<https://github.com/attractivechaos/klib>), distributed under the MIT License;
the license text is at the top of the file.

## sse2neon

`csrc/dp/sse2neon/emmintrin.h` is the sse2neon SSE-to-NEON translation header
(<https://github.com/DLTcollab/sse2neon>), used on ARM builds only and
distributed under the MIT License; the license text is in the file.

## zlib

zlib is not part of this repository; builds link the zlib found on the build
host. The Linux binary release links zlib 1.3.2, built from the release tarball
at <https://github.com/madler/zlib> (pinned by checksum), statically.

zlib is distributed under the zlib License:

> Copyright (C) 1995-2026 Jean-loup Gailly and Mark Adler
>
> This software is provided 'as-is', without any express or implied
> warranty. In no event will the authors be held liable for any damages
> arising from the use of this software.
>
> Permission is granted to anyone to use this software for any purpose,
> including commercial applications, and to alter it and redistribute it
> freely, subject to the following restrictions:
>
> 1. The origin of this software must not be misrepresented; you must not
>    claim that you wrote the original software. If you use this software
>    in a product, an acknowledgment in the product documentation would be
>    appreciated but is not required.
>
> 2. Altered source versions must be plainly marked as such, and must not be
>    misrepresented as being the original software.
>
> 3. This notice may not be removed or altered from any source distribution.
