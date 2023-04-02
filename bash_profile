
# this makes your prompt much cooler
user_color=36
host_color=32
path_color="33;1"
export PS1="\[\033["$user_color"m\]\u\[\033[m\]@\[\033["$host_color"m\]\h:\[\033["$path_color"m\]\w\[\033[m\]\$ "

# this colors your files and folders differently
export LS_OPTIONS='--color=auto'
eval "$(dircolors -b)"



alias ls='ls $LS_OPTIONS'
