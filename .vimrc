" VIM SETTING {{{
set nu
set hlsearch
set incsearch
set showcmd
set foldmethod=marker
syntax on
filetype plugin indent on
set path+=~/,**,*,;
packadd! matchit
" }}}

" VIM MAPPING {{{
let mapleader=','
let maplocalleader=';'

" NO MORE ESC KEY
inoremap jk <esc>
inoremap <esc> <nop>
nnoremap <leader>sv :w!<cr> :source ~/.vimrc<cr>
nnoremap <leader>ev :w!<cr> :vsplit ~/.vimrc<cr>
nnoremap <localleader>py :w!<cr> :!clear<cr> :!python %<cr>

" Bracket around word
nnoremap <leader>wb mai{<esc>ea}<esc>`a
nnoremap <leader>wp mai(<esc>ea)<esc>`a
" }}}

" VIM PLUGIN {{{
execute pathogen#infect()
"  }}}

" VIM PLUGIN SETTING {{{
" make YCM compatible with UltiSnips (using supertab)
let g:ycm_key_list_select_completion = ['<C-n>', '<Down>']
let g:ycm_key_list_previous_completion = ['<C-p>', '<Up>']
let g:SuperTabDefaultCompletionType = '<C-n>'

" better key bindings for UltiSnipsExpandTrigger
let g:UltiSnipsExpandTrigger = "<tab>"
let g:UltiSnipsJumpForwardTrigger = "<tab>"
let g:UltiSnipsJumpBackwardTrigger = "<s-tab>"

"  }}}
