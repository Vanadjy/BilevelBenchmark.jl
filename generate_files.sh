module load julia/1.11

#purposes=("Known_Referee" "Anonymous_Referee")
purposes=("Compare_Lambda")
for purpose in "${purposes[@]}"; do
   if [[ "$purpose" == "Known_Referee" ]]; then
      refs=("Intern_EndPoint" "Intern_Complete" "Intern_Reverse")
   elif [[ "$purpose" == "Anonymous_Referee" ]]; then
      refs=("Extern_EndPoint" "Extern_Complete" "Extern_Reverse")
   elif [[ "$purpose" == "Compare_Lambda" ]]; then
      refs=()
      effort_choices=("Agregate" "UL" "LL")
      for effort_choice in "${effort_choices[@]}"; do
         julia --project=. numerics/ProfilesHub.jl "dec" "Extern_EndPoint" $purpose $effort_choice "false";
      done
   else
      echo "Invalid purpose. Please set the purpose variable to either 'Known_Referee', 'Anonymous_Referee', or 'Compare_Lambda'."
   fi
   for ref in "${refs[@]}"; do
      julia --project=. numerics/ProfilesHub.jl "dec" $ref $purpose "Agregate" "true";
   done
done