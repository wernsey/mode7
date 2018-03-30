# Mode 7

[Mode7-style][wiki] graphics demo for SDL

## Compiling 

It has only been tested under Windows with the MinGW tools.
The `Makefile` builds it for Windows with SDL 2. 
`Makefile.gdi` builds it for Windows using GDI - that is, without any third-party dependencies.

## References 

The math is all derived from Jasper _"cearn"_ Vijn's [TONC book][tonc] on Gameboy Advanced programming;
Chapter 21 has the description of all the relevant algorithms.

* The demo uses graphics from the [Hard Vacuum tile set][vacuum]
* The sprite comes from [charas-project.net][charas]

The demo uses my own [bitmap][] module.

[wiki]: https://en.wikipedia.org/wiki/Mode_7
[tonc]: http://www.coranac.com/tonc/text/mode7ex.htm
[vacuum]: http://www.lostgarden.com/2005/03/game-post-mortem-hard-vacuum.html
[charas]: http://charas-project.net/
[bitmap]: https://github.com/wernsey/bitmap