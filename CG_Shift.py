# This is just a test
#
import glob
import os
#import matplotlib.pyplot as plt
#import numpy as np

# Find all files with the extension '.inp'
files = glob.glob("*.inp")

# Read and execute each input file
for file in files:
    name = os.path.splitext(os.path.basename(file))[0]
    with open(file, 'r') as input:
        code = input.read()
        exec(code)

    # Define CG functions to be compared against target CGs:
    def XX(W1, W2, W3):
        return (W0*X0+W1*X1+W2*X2+W3*X3)/(W0+W1+W2+W3)
    def YY(W1, W2, W3):
        return (W0*Y0+W1*Y1+W2*Y2+W3*Y3)/(W0+W1+W2+W3)
    def ZZ(W1, W2, W3):
        return (W0*Z0+W1*Z1+W2*Z2+W3*Z3)/(W0+W1+W2+W3)

    # Cycle through adjustments m times, n adjustments each cycle:
    n = 96
    m = 48
    W1A = W1B = W10 = W1
    W2A = W2B = W20 = W2
    W3A = W3B = W30 = W3
    #iterations = np.array([])
    #W1_points = np.array([])
    #W2_points = np.array([])
    #W3_points = np.array([])
    #k = 0
    for i in range(n):
        for j in range(m):
            W1A = W1 + W10/(2**(j+i))
            W1B = W1 - W10/(2**(j+i))
            if (abs(X-XX(W1A, W2, W3)) < abs(X-XX(W1B, W2, W3))):
                W1 = W1A
            else:
                W1 = W1B
            #k = k+1
            #iterations = np.append(iterations, k)
            #W1_points = np.append(W1_points, W1)
        for j in range(m):
            W2A = W2 + W20/(2**(j+i))
            W2B = W2 - W20/(2**(j+i))
            if (abs(Y-YY(W1, W2A, W3)) < abs(Y-YY(W1, W2B, W3))):
                W2 = W2A
            else:
                W2 = W2B
            #W2_points = np.append(W2_points, W2)
        for j in range(m):
            W3A = W3 + W30/(2**(j+i))
            W3B = W3 - W30/(2**(j+i))
            if (abs(Z-ZZ(W1, W2, W3A)) < abs(Z-ZZ(W1, W2, W3B))):
                W3 = W3A
            else:
                W3 = W3B
            #W3_points = np.append(W3_points, W3)

    # Weight scaling factors
    fW = W/(W0+W1+W2+W3)
    f1 = fW*W1/W10
    f2 = fW*W2/W20
    f3 = fW*W3/W30

    # Print into terminal
    print()
    print("From '{}.inp'".format(name))
    print('------------------------------------')
    print('Group weight factors:')
    print()
    print('SELFWEIGHT Y {:.3f}'.format(round(-fW,3)),Name_0)
    print('SELFWEIGHT Y {:.3f}'.format(round(-f1,3)),Name_X)
    print('SELFWEIGHT Y {:.3f}'.format(round(-f2,3)),Name_Y)
    print('SELFWEIGHT Y {:.3f}'.format(round(-f3,3)),Name_Z)
    print()
    print('Acceleration factors:')
    print()
    print('SELFWEIGHT X {:.3f}'.format(round(fW*0.05,3)),Name_0)
    print('SELFWEIGHT X {:.3f}'.format(round(f1*0.05,3)),Name_X)
    print('SELFWEIGHT X {:.3f}'.format(round(f2*0.05,3)),Name_Y)
    print('SELFWEIGHT X {:.3f}'.format(round(f3*0.05,3)),Name_Z)
    print('*SELFWEIGHT X {:.3f}'.format(round(0.05)),'_...')
    print()
    print('SELFWEIGHT Z {:.3f}'.format(round(fW*0.05,3)),Name_0)
    print('SELFWEIGHT Z {:.3f}'.format(round(f1*0.05,3)),Name_X)
    print('SELFWEIGHT Z {:.3f}'.format(round(f2*0.05,3)),Name_Y)
    print('SELFWEIGHT Z {:.3f}'.format(round(f3*0.05,3)),Name_Z)
    print('*SELFWEIGHT Z {:.3f}'.format(round(0.05)),'_...')
    print()

    # Print into output file
    filename = name+'.out'
    if (Outputfile == True):
        with open(filename, 'w') as output:
            output.write('\n')
            output.write("From '{}.inp'\n".format(name))
            output.write('------------------------------------'+'\n')
            output.write('Group weight factors:'+'\n')
            output.write('\n')
            output.write('SELFWEIGHT Y {:.3f} '.format(-fW)+Name_0+'\n')
            output.write('SELFWEIGHT Y {:.3f} '.format(-f1)+Name_X+'\n')
            output.write('SELFWEIGHT Y {:.3f} '.format(-f2)+Name_Y+'\n')
            output.write('SELFWEIGHT Y {:.3f} '.format(-f3)+Name_Z+'\n')
            output.write('\n')
            output.write('Acceleration factors:'+'\n')
            output.write('\n')
            output.write('SELFWEIGHT X {:.3f} '.format(fW*0.05)+Name_0+'\n')
            output.write('SELFWEIGHT X {:.3f} '.format(f1*0.05)+Name_X+'\n')
            output.write('SELFWEIGHT X {:.3f} '.format(f2*0.05)+Name_Y+'\n')
            output.write('SELFWEIGHT X {:.3f} '.format(f3*0.05)+Name_Z+'\n')
            output.write('SELFWEIGHT X {:.3f} '.format(0.05)+'_...'+'\n')
            output.write('\n')
            output.write('SELFWEIGHT Z {:.3f} '.format(fW*0.05)+Name_0+'\n')
            output.write('SELFWEIGHT Z {:.3f} '.format(f1*0.05)+Name_X+'\n')
            output.write('SELFWEIGHT Z {:.3f} '.format(f2*0.05)+Name_Y+'\n')
            output.write('SELFWEIGHT Z {:.3f} '.format(0.05)+'_...'+'\n')

    # Print plots
    #plt.plot(iterations, W1_points)
    #plt.plot(iterations, W2_points)
    #plt.plot(iterations, W3_points)
    #plt.show()
