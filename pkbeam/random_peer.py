#%%
import numpy as np 
import random 
bravo = ["Carlos","Dave","Sebbe","Guillaume","Yuning"]
paper1 = bravo
paper2 = bravo
random.shuffle(paper1)
random.shuffle(paper2)
review = {}
for name in bravo:
    review[name]=[] 
for person in bravo:
    for report1 in paper1:
        if report1 not in review[person] and person != report1 and report1 not in review[person] :
            review[person].append(report1)
        for report2 in paper2:
            if person !=report2 and report1 != report2 and report2 not in review[person]:
                review[person].append(report2)

final = {}
for name in bravo:
    final[name]=[] 

discard = []
for name in bravo:
    print("Assigning for {}:".format(name))
    for i in range(2):
        
        peer1 = random.choice(review[name])
        # print(peer1)
        print("\tThe paper {} is from {}".format(i+1,peer1))
        final[name].append(peer1)

        if peer1 in review[name]:
            review[name].remove(peer1)
        
        discard.append(peer1)

    for ele in discard:
            if discard.count(ele) ==2 :
                for name in bravo:
                    if ele in review[name]:
                        review[name].remove(ele)

print(final)
import pandas as pd
assign = pd.DataFrame(final)
