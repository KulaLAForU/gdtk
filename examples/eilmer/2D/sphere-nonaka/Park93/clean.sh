function delete_items() {
    for item in "$@"; do
        if [[ -f "$item" ]]; then
            rm "$item"  # Deletes a file
        elif [[ -d "$item" ]]; then
            rm -rf "$item"  # Deletes a folder and its contents
        fi
    done
}
printf "Everything's been cleaned up...\n"
# Example usage
delete_items "./config" "./flow" "./grid" "./hist" "./loads" "./plot" "./solid" *.gas *.chem  *.pdf *.dat *.exch *.csv "./LOGFILE"
