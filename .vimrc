" Vim Settings ---------------------{{{
set nocompatible              " be iMproved, required
set backspace=indent,eol,start
set foldmethod=marker
set number
set ruler
set showcmd
set incsearch
set hlsearch
set wildmenu
set wildmode=longest:full,full
set tabstop=3
set shiftwidth=3
set expandtab
set smartindent
colorscheme murphy
filetype off                  " required
" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
set rtp+=~/.vim/plugin/
set background=dark
syntax enable
set encoding=utf-8

let g:netrw_use_noswf= 0

let maplocalleader = ','
let mapleader = ';'

"}}}

" Plugins --------------------{{{
execute pathogen#infect()
syntax on
filetype plugin indent on

" matchit
packadd! matchit
"}}}

" Plugin Settings --------------------{{{

" Termcolor setting {{{
let g:solarized_termcolors=256

if has('gui_running')
	set background=light
else
	set background=dark
endif
" Color setting }}}

" Powerline setup {{{
set guifont=DejaVu\ Sans\ Mono\ for\ Powerline\ 9
set laststatus=2
" }}}

" dbext MySQL default profile {{{
let g:dbext_default_profile_temp_mysql = 'type=MYSQL:user=root:passwd=CkMj1527!%@&:dbname=temporary:host=localhost:extra=-t'
let g:dbext_default_profile_employee = 'type=MYSQL:user=root:passwd=CkMj1527!%@&:dbname=employees:host=localhost:extra=-t'

let g:dbext_default_profile_mySQL = 'type=MYSQL:user=root:passwd=CkMj1527!%@&:dbname=temporary'  
let g:dbext_default_window_use_horiz = 0  " Use vertical split
let g:dbext_default_window_width = 70

let g:C_UseTool_cmake ='yes'
let g:C_UseTool_doxygen = 'yes'
" }}}

" Nvim-R {{{
" always vertical split
let R_rconsole_width = 57
let R_min_editor_width = 18
let R_in_buffer=1
let R_pdfviewer = 'zathura'
let R_complete = 2 
" Nvim-R Send line
vmap <Space> <Plug>RDSendSelection
nmap <Space> <Plug>RDSendLine
" vimtex mupdf viewing pdf through synctex
let g:vimtex_view_method = 'mupdf'
" vimtex start server
"if empty(v:servername) && exists('*remote_startserver')
    "call remote_startserver('VIM')
"endif
" }}}

" python3 syntax checkin {{{
let g:pymode_python = 'python3'

" Run python in vim
let g:pymode_run = 1
let g:pymode_run_bind = '<leader>r'
" }}}

" YCM {{{
let g:ycm_global_ycm_extra_conf = '~/.vim/bundle/YouCompleteMe/third_party/ycmd/.ycm_extra_conf.py'
let g:ycm_python_binary_path = 'python3.6'
let g:C_CplusCFlags = '-Wall -g -o0 -c -std=c++1z -I/home/ck1/Downloads/boost_1_67_0 -I../include/armadillo -DARMA_DONT_USE_WRAPPER -lopenblas -llapack'
let g:C_CplusLFlags = '-Wall -g -o0 -I/home/ck1/Downloads/boost_1_67_0/libs -I../include/armadillo -DARMA_DONT_USE_WRAPPER -lopenblas -llapack '
let g:C_CplusLibs = '-llapack -lopenblas'
" make YCM compatible with UltiSnips (using supertab)
let g:ycm_key_list_select_completion = ['<C-n>', '<Down>']
let g:ycm_key_list_previous_completion = ['<C-p>', '<Up>']
let g:SuperTabDefaultCompletionType = '<C-n>'

" better key bindings for UltiSnipsExpandTrigger
let g:UltiSnipsExpandTrigger = "<tab>"
let g:UltiSnipsJumpForwardTrigger = "<tab>"
let g:UltiSnipsJumpBackwardTrigger = "<s-tab>"
" }}}

" Vim-GDB {{{
let g:ConqueTerm_Color = 2
let g:ConqueTerm_CloseOnEnd = 1
let g:ConqueTerm_StartMessages = 0
" }}}

" AsyncRun quickfix {{{
let g:asyncrun_open = 8
" }}}

" FZF {{{
set rtp+=~/.fzf
" }}}

"}}}

" Vim Autocommands --------------------{{{
augroup vimrc_autocmds
	    autocmd!
        " highlight characters past column 120
        autocmd FileType python highlight Excess ctermbg=DarkGrey guibg=Black
        autocmd FileType python match Excess /\%120v.*/
        autocmd FileType python set nowrap
augroup END

augroup runcodegroup
	    autocmd!
        autocmd FileType python nnoremap <buffer> <D-r> :Pyrun<CR><C-W>jG<C-W>k
augroup END

" File Types you want to use suggestions
autocmd FileType php setlocal completefunc=MySQLCompleteFunction
autocmd FileType javascript setlocal completefunc=MySQLCompleteFunction
autocmd FileType sql setlocal completefunc=MySQLCompleteFunction
" cppman
command! -nargs=+ Cppman silent! call system("tmux split-window cppman " . expand(<q-args>))
autocmd FileType cpp nnoremap <silent><buffer> K <Esc>:Cppman <cword><CR>

augroup template
    autocmd!
    autocmd BufNewFile *.tex 0r ~/.vim/templates/skeleton.tex |ks|call FirstCreate()|'s
    autocmd BufNewFile *.py 0r ~/.vim/templates/skeleton.py |ks|call FirstCreate()|'s
    autocmd BufNewFile *.R,*.r 0r ~/.vim/templates/skeleton.R |ks|call FirstCreate()|'s
    autocmd BufNewFile *.Rnw,*.rnw 0r ~/.vim/templates/skeleton.Rnw |ks|call FirstCreate()|'s
    fun FirstCreate()
        if line("$") > 10
            let l = 10
        else
            let l = line("$")
        endif
        exe "1," . l . "g/Filename      : /s/Filename      : .*/Filename      : " . expand("%:p:t")
        exe "1," . l . "g/Date created  : /s/Date created  : .*/Date created  : " . strftime("%c")
        exe "1," . l . "g/Created by    : /s/Created by    : .*/Created by    : " . $USER
    endfun

    autocmd BufWritePre,FileWritePre *.py,*.R,*.r,*.Rnw,*.rnw ks|call LastMod()|'s
    fun LastMod()
        if line("$") > 10
            let l = 10
        else
            let l = line("$")
        endif
        exe "1," . l . "g/Last modified : /s/Last modified : .*/Last modified : " . strftime("%c")
        exe "1," . l . "g/Modified by   : /s/Modified by   : .*/Modified by   : " . $USER
    endfun
augroup END


"}}}

" Vim Custom Functions --------------------{{{
" ex command for toggling hex mode - define mapping if desired
if !exists("*ToggleBinmode")
    command -bar HexBinmode call ToggleBinmode()

    function ToggleBinmode()
        " Bin mode should be considered a read-only operation
        " save values for modified and read-only for restoration later,
        " and clear the read-only flag for now
        let l:modified=&mod
        let l:oldreadonly=&readonly
        let &readonly=0
        let l:oldmodifiable=&modifiable
        let &modifiable=1
        if !exists("b:editHex") || !b:editHex
            " save old options
            let b:oldft=&ft
            let b:oldbin=&bin
            " set new options
            setlocal binary
            silent :e
            let &ft="xxd"
            let b:editHex=1
            %!xxd -b
        else
            let &ft=b:oldft
            if !b:oldbin
                setlocal binary
            endif
            let b:editHex=0
            %!xxd -r
        endif
        let &mod=l:modified
        let &readonly=l:oldreadonly
        let &modifiable=l:oldmodifiable
    endfunction
endif

if !exists("*ToggleHex")
    command -bar Hexmode call ToggleHex()
    function ToggleHex()
        " hex mode should be considered a read-only operation
        " save values for modified and read-only for restoration later,
        " and clear the read-only flag for now
        let l:modified=&mod
        let l:oldreadonly=&readonly
        let &readonly=0
        let l:oldmodifiable=&modifiable
        let &modifiable=1
        if !exists("b:editHex") || !b:editHex
            " set new options
            setlocal binary " make sure it overrides any textwidth, etc.
            silent :e " this will reload the file without trickeries 
            "(DOS line endings will be shown entirely )
            let &ft="xxd"
            " set status
            let b:editHex=1
            " switch to hex editor
            %!xxd
        else
            " restore old options
            let &ft=b:oldft
            if !b:oldbin
                setlocal nobinary
            endif
            " set status
            let b:editHex=0
            " return to normal editing
            %!xxd -r
        endif
        " restore values for modified and read only state
        let &mod=l:modified
        let &readonly=l:oldreadonly
        let &modifiable=l:oldmodifiable
    endfunction
endif

if !exists("*DebugSession")
    function DebugSession()
        silent make -o vimgdb -gcflags "-N -l"
        redraw!
        if (filereadable("vimgdb"))
            ConqueGdb vimgdb
        else
            echom "Couldn't find debug file"
        endif
    endfunction
endif

if !exists("*DebugSessionCleanup")
    function DebugSessionCleanup(term)
        if (filereadable("vimgdb"))
            let ds=delete("vimgdb")
        endif
    endfunction
    call conque_term#register_function("after_close", "DebugSessionCleanup")
endif
nmap <leader>D :call DebugSession()<CR>;
" foldcolun toggle
function! FoldColumnToggle()
    if &foldcolumn
        setlocal foldcolumn=0
    else
        setlocal foldcolumn=4
    endif
endfunction

" quickfix togle
let g:quickfix_is_open = 0
function! Quickfixtoggle()
    if g:quickfix_is_open
        cclose
        let g:quickfix_is_open = 0
        exe g:quickfix_return_to_window . "wincmd w"
    else
        let g:quickfix_return_to_window = winnr()
        copen
        let g:quickfix_is_open = 1
    endif
endfunction

" AsyncRun Compile {{{
function! s:compile_and_run()
    exec 'w'
    if &filetype == 'c'
        exec "AsyncRun! gcc % -o %<; time ./%<"
    elseif &filetype == 'cpp'
        exec "AsyncRun! g++ -std=c++11 % -o %<; time ./%<"
    elseif &filetype == 'java'
        exec "AsyncRun! javac %; time java %<"
    elseif &filetype == 'sh'
        exec "AsyncRun! time bash %"
    elseif &filetype == 'python'
        exec "AsyncRun! time python %"
    endif
endfunction
" }}}

" Rtags creation {{{
command! Rtags :AsyncRun! Rscript -e 'library(here); rtags(path=paste0(here()), recursive=TRUE, ofile="RTAGS")' -e 'etags2ctags("RTAGS", "~/rtagspkg")' -e 'unlink("RTAGS")' -e 'rtags(path=paste0("~/R/x86_64-pc-linux-gnu-library/3.6"), recursive=TRUE, ofile="RTAGS")' -e 'etags2ctags("RTAGS", "~/rtags")' -e 'unlink("RTAGS")'

" }}}

" }}}

" Mappings --------------------{{{

" Custom Mappings {{{
inoremap <leader>stf <esc>:w!<cr>
nnoremap <leader>stf :w!<cr>
inoremap <c-u> <esc>bviwUA
nnoremap <leader>ev :vsplit $MYVIMRC<cr>G
nnoremap <leader>sv :w!<cr> :so $MYVIMRC<cr>
nnoremap -dd ddo
nnoremap -dw bdiwa
nnoremap <leader>" viw<esc>a"<esc>bi"<esc>lel
nnoremap <leader>' viw<esc>a'<esc>bi'<esc>lel
inoremap jk <esc>
inoremap <esc> <nop>
onoremap ih :<c-u>execute "normal! ?\\(^==\\+\\\|^--\\+\\)$\r:nohlsearch\rkvg_"<cr>
nmap <C-N><C-N> :set invnumber<CR>
nnoremap <LocalLeader>gf <Plug>RSyncFor
nnoremap <leader>cd :cd %:p:h<cr>
" grep
":nnoremap <leader>g :silent exe "grep! -R " . shellescape(expand("<cWORD>")) . " ."<cr>:copen<cr>
" Number
nnoremap <leader>sn :setlocal number!<cr>
" foldcolumn toggle
nnoremap <leader>sf :call FoldColumnToggle()<cr>
" quickfix toggle
nnoremap <leader>qt :call Quickfixtoggle()<cr>
nnoremap <F4> :call <SID>compile_and_run()<CR>
" }}}

" Split screen navigation {{{
nnoremap <C-J> <C-W><C-J>
nnoremap <C-K> <C-W><C-K>
nnoremap <C-L> <C-W><C-L>
nnoremap <C-H> <C-W><C-H>
" }}}

" vim HexMode {{{
nnoremap <localleader>hm :Hexmode<CR>
inoremap <localleader>hm <Esc>:Hexmode<CR>
vnoremap <localleader>hm :<C-U>Hexmode<CR>

nnoremap <C-B> :HexBinmode<CR>
inoremap <C-B> <Esc>:HexBinmode<CR>
vnoremap <C-B> :<C-U>HexBinmode<CR>
" }}}

" " For vim-fugitive {{{
nnoremap <space>ga :Git add %<CR><CR>
nnoremap <space>gs :Gstatus<CR>
nnoremap <space>gc :Gcommit -v -q<CR>
nnoremap <space>gt :Gcommit -v -q %:p<CR>
nnoremap <space>gd :Gdiff<CR>
nnoremap <space>ge :Gedit<CR>
nnoremap <space>gr :Gread<CR>
nnoremap <space>gw :Gwrite<CR><CR>
nnoremap <space>gl :silent! Glog<CR>:bot copen<CR>
nnoremap <space>gp :Ggrep<Space>
nnoremap <space>gm :Gmove<Space>
nnoremap <space>gb :Git branch<Space>
nnoremap <space>gbl :Gblame<CR>
nnoremap <space>go :Git checkout<Space>
nnoremap <space>gps :Dispatch! git push origin master<CR>
nnoremap <space>gpsb :Dispatch! git push origin 
nnoremap <space>gpl :Dispatch! git pull origin master<CR> " vim tabs
nnoremap <C-Left> :tabprevious<CR>
nnoremap <C-Right> :tabnext<CR>
nnoremap <silent> <A-Left> :execute 'silent! tabmove ' . (tabpagenr()-2)<CR>
nnoremap <silent> <A-Right> :execute 'silent! tabmove ' . (tabpagenr()+1)<CR>
" }}}

" match trailing white spaces as error {{{
nnoremap <leader>dw :exe "match Error /" . '\v( )+$/'<CR>
nnoremap <leader>rv /\v
" }}}


" vim fzf {{{
nnoremap <leader>ff :Files<cr> " fuzzy find files in working directory (where vim launched)"
nnoremap <leader>f/ :BLines<cr> "fuzzy find lines in the current file"
nnoremap <leader>fb :Buffers<cr> "fuzzy find an open buffer
nnoremap <leader>fr :Rg<cr> "funnzy find text in the working directory"
nnoremap <leader>fc :Commands<cr> " fuzzy find Vim commands (like Ctr-Shift-P in sublime/atom/vsc)"
" }}}

" }}}

" Abbreviations --------------------{{{
:iabbrev adn and
:iabbrev waht what
:iabbrev tehn then
:iabbrev <expr> dts strftime("%c")
:iab @@ chong.kim@ucdenver.edu
"}}}

" Tags {{{
augroup Tags
   autocmd!
   autocmd FileType r set tags+=~/.cache/Nvim-R/Rtags,~/.cache/Nvim-R/RsrcTags,~/R/x86_64-pc-linux-gnu-library/3.6/tags,~/R/x86_64-pc-linux-gnu-library/3.6/Rtags,*/rtags,**/rtags,~/rtags,~/rtagspkg
   autocmd FileType rnoweb set tags+=~/.cache/Nvim-R/Rtags,~/.cache/Nvim-R/RsrcTags,~/R/x86_64-pc-linux-gnu-library/3.6/tags,~/R/x86_64-pc-linux-gnu-library/3.6/Rtags,*/rtags,**/rtags,~/rtags,~/rtagspkg

augroup END

let g:tagbar_type_r = {
    \ 'ctagstype' : 'r',
    \ 'kinds'     : [
        \ 'f:Functions',
        \ 'g:GlobalVariables',
        \ 'v:FunctionVariables',
    \ ]
    \ }


" }}}
