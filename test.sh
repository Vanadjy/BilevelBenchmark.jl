purposes=("Anonymous_Referee" "Known_Referee" "Compare_Lambda")
efforts=("UL" "Agregate")
for purpose in "${purposes[@]}"; do
    echo "Running ProfilesHub for purpose: $purpose"
    if [ $purpose = "Anonymous_Referee" ]; then
        echo "Running with external referees"
        refs=("Extern_EndPoint" "Extern_Complete" "Extern_Reverse"); 
    else
        echo "Running with internal referees"
        refs=("Intern_EndPoint" "Intern_Complete" "Intern_Reverse"); 
    fi; 
    for ref in "${refs[@]}"; do
        for effort in "${efforts[@]}"; do
            julia --project=. numerics/ProfilesHub.jl "dec" $ref $purpose $effort;
        done;
    done;
done