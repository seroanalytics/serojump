#!/bin/bash

BRANCH_NAME="main"
CHUNK_SIZE=10

# Get the list of commits to push
commits=$(git rev-list --reverse HEAD)

# Convert commits to an array
commit_array=($commits)

# Push commits in chunks
for ((i=0; i<${#commit_array[@]}; i+=CHUNK_SIZE)); do
    # Determine the end of the chunk
    end=$((i+CHUNK_SIZE-1))
    if [ $end -ge ${#commit_array[@]} ]; then
        end=$((${#commit_array[@]}-1))
    fi

    # Get the range of commits to push
    start_commit=${commit_array[$i]}
    end_commit=${commit_array[$end]}

    # Push the range of commits
    git push origin $end_commit:refs/heads/$BRANCH_NAME --force-with-lease

    echo "Pushed commits from $start_commit to $end_commit"
done