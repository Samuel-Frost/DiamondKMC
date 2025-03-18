### TODO:
- Generate cells with custom distributions for different defects e.g. evenly distributed Nitrogen, I linked with V
- Make extxyz output continuous instead of being dumped at the end
- When vacancies and interstitials combine, they can just both have the same pos, don't write to extxyz
- Dissociation needs to be implemented -> how to pick the direction they move in
- For vacancies, the formation of di/tri/quad is directional, needs to be implemented
- Implement strain fields, make it so there is a bias vector field, e.g. can have a divacancy have a negative bias anywhere but inplane 
- Now, rates need to be picked from basically what is possible in every direction for a defect
- Store defect types not in generic defect class but as inherited classes, then vacancy chain can also be its own class, easier to enforce rules
- ^^ no stupid combine function

1. initialise every defect to have pairs = [their index]
2. when two defects combine, set one defect to have their pair list be *the same* as the other (a pointer), and then add the other index to it 
3. now both defect share literally the same place in memory where their pair list is
4. when another defect wants to join, make its pair list equal the others, then add it to the list
5. when one dissociates, remove it from the single pair  list, and then make its pair list its index again.
