#Coin flip streaks
import random
numberOfStreaks = 0
for experimentNumber in range(10000):
    #Code that creates a list of 100 'heads' or 'tails' values.
    results = ''
    for i in range (100):
        flip = random.randint(0,1)
        results += str(flip)
   
    #Code that checks if there is a streak of 6 heads or tails in a row.
    if '111111' in results:
        numberOfStreaks += 1
        results = ''
    elif '000000' in results:
        numberOfStreaks += 1
        results = ''
    else:
        numberOfStreaks += 0
        results = ''
        
print('Number of experiments with 9-streak: ' + str(numberOfStreaks))
print('Experiments run: ' + str(experimentNumber+1))
print('Chance of streak: %s%%' % (numberOfStreaks / 100))
