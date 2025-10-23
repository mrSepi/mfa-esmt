import heapq
import itertools
import math
import random
import resource
import statistics
import subprocess
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

def emst(a):
    # restituisce matrice adiacenza
    matrix = []
    for i in range(len(a)):
        matrix.append([0 for k in range(len(a))])
    index = {}
    for i in range(len(a)):
        index[a[i]] = i

    REMOVED = 'removed'
    counter = itertools.count()
    res = []
    info = {}
    queue = []
    for p in a:
        entry = [math.inf, next(counter), p, None]
        info[p] = entry
        queue.append(entry)

    while info:
        d, _, x, v = heapq.heappop(queue)
        if x is not REMOVED:
            del info[x]
            if not math.isinf(d):
                res.append((x,v))
                matrix[index[x]][index[v]] = 1
                matrix[index[v]][index[x]] = 1
            for _, y in info.items():
                dist = math.dist(x, y[2])
                if dist < y[0]:
                    changed_priority = [dist, next(counter), y[2], x]
                    info[y[2]] = changed_priority
                    y[2] = REMOVED
                    heapq.heappush(queue, changed_priority)
    return matrix
def draw(points, color):
      plt.scatter([p[0] for p in points], [p[1] for p in points], c=color)
def drawBig(points, color, times):
    plt.scatter([p[0] for p in points], [p[1] for p in points], c=color, s=[20*4**times for i in range(len(points))])
def drawSegment(p1, p2) :
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], c='darkgray')
def drawCircle(r):
    xs = []
    ys = []
    angoli = np.arange(0,2*math.pi+0.1,0.1)
    for alpha in angoli:
        xs.append(math.cos(alpha)*r)
        ys.append(math.sin(alpha)*r)
    plt.plot(xs, ys, c='lightgray', linestyle='dashed')
def drawAnnotate(points, color, left=False):
      plt.axis('equal')
      plt.grid(True)
      plt.scatter([p[0] for p in points], [p[1] for p in points], c=color)
      coord_label = (4,4) if not left else (-4,4)
      for i in range(len(points)):
        plt.annotate("{}".format(i), points[i], textcoords="offset points", xytext=coord_label, c=color)
def drawPointWithDegree(points, color):
      plt.axis('equal')
      plt.grid(True)
      plt.scatter([p[0] for p in points], [p[1] for p in points], c=color)
      for i in range(len(points)):
        deg = 0
        for t in range(n):
            if x[i][t] == 1:
                deg += 1
        for s in range(n-2):
            if y[i][s] == 1:
                deg += 1
        plt.annotate("{}".format(deg), points[i], textcoords="offset points", xytext=(4,4), c='black')
def drawHeu(a, sp, radius=0):
    draw(sp, 'red')
    draw(a, 'black')
    points = a + sp
    matrix = emst(points)
    for t in range(len(points)):
      for s in range(len(points)):
        if matrix[s][t] == 1:
          drawSegment(points[t], points[s])
    if radius > 0:
        drawCircle(radius)

    fig = plt.gcf()
    fig.set_size_inches(8,6)
    #plt.suptitle("vs optimal: +%.3f%%" % (((len_heu_best/len_emst(points))-1)*100))
    plt.savefig('/tmp/heu.png', bbox_inches='tight')
def drawOpt(a, sp, optimalsp):
    axes = plt.subplot(1,2,1)
    axes.clear()
    plt.axis('equal')
    plt.grid(True)
    draw(sp, 'red')
    draw(a, 'black')
    points = a + sp
    matrix = emst(points)
    for t in range(len(points)):
      for s in range(len(points)):
        if matrix[s][t] == 1:
          drawSegment(points[t], points[s])

    axes = plt.subplot(1,2,2)
    axes.clear()
    plt.axis('equal')
    plt.grid(True)
    drawAnnotate(a, 'black')
    drawAnnotate(optimalsp, 'red')
    points = a + optimalsp
    matrix = emst(points)
    for t in range(len(points)):
      for s in range(len(points)):
        if matrix[s][t] == 1:
          drawSegment(points[t], points[s])

    fig = plt.gcf()
    fig.set_size_inches(16,6)
    # plt.suptitle("vs optimal: +%.3f%%" % surplus)
    plt.savefig('/tmp/opt.png', bbox_inches='tight')

def printMatrixSideBySide(x, y) :
    for i in range(len(x)):
        for j in range(len(x[i])):
            print('%.3f' % x[i][j], end=" ")
        print(" |", end = " ")
        for j in range(len(y[i])):
            print('%.3f' % y[i][j], end=" ")
        print()
def drawSegmentColoured(p1, p2, color, alpha=1) :
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], c=color, alpha=alpha)
def drawSegmentColouredBig(p1, p2, color, width, alpha=1) :
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], c=color, alpha=alpha, lw=width)

def drawHeu12(a, sp, x, y):
    n = len(a)

    axes = plt.subplot(1,2,1)
    axes.clear()
    plt.axis('equal')
    plt.grid(True)
    drawAnnotate(sp, 'red', True)
    drawAnnotate(a, 'black')
    for t in range(n):
        max_val = 0
        index = 0
        for s in range(n-2):
          if x[s][t] > max_val:
            max_val = x[s][t]
            index = s
        drawSegmentColoured(a[t], sp[index], 'black')
    for s in range(n-2):
      for r in range(n-2):
        if y[s][r] >= 1:
          drawSegmentColoured(sp[r], sp[s], 'red')

    axes = plt.subplot(1,2,2)
    axes.clear()
    plt.axis('equal')
    plt.grid(True)
    drawAnnotate(sp, 'red')
    drawAnnotate(a, 'black')
    points = a + sp
    matrix = emst(points)
    for t in range(len(points)):
      for s in range(len(points)):
        if matrix[s][t] == 1:
          drawSegment(points[t], points[s])
    
    fig = plt.gcf()
    fig.set_size_inches(16,6)
    plt.savefig('/tmp/heu12.png', bbox_inches='tight')
def normalizeMatrix(x, n):
    for t in range(n):
        max_val = 0
        index = 0
        for s in range(n-2):
          if x[s][t] > max_val:
            x[index][t] = 0
            max_val = x[s][t]
            index = s
        if max_val != 1:
            print("x[][] = %.3f\n" % max_val, file=sys.stderr)
        x[index][t] = 1
def get_steiner_points_and_matrix(f_name):
    x  = []
    y  = []
    sp = []

    file = open(f_name, "r")
    readingMatrix  = True
    line = file.readline()
    while line:
        if readingMatrix:
            print(line)
        if line == "\n":
            readingMatrix = False
        elif readingMatrix:
            xs, ys = line.split("|")
            x.append([float(c) for c in xs.strip().split(" ")])
            y.append([float(c) for c in ys.strip().split(" ")])
        else:
            a1, a2 = line.split(",")
            sp.append((float(a1),float(a2)))
        line = file.readline()
    file.close()
    printMatrixSideBySide(x,y)
    return x, y, sp
def draw_instance2(finput_name, exe_name, fsol_name="", fsol_has_lenghts=True, radius=0):
    a = get_terminals(finput_name)

    input_file  = open(finput_name, "r")
    output_file = open("Output.txt", "w")
    subprocess.run(["./" + exe_name], stdin=input_file, stdout=output_file)
    input_file.close()
    output_file.close()
    x, y, sp = get_steiner_points_and_matrix("Output.txt")
    normalizeMatrix(x, len(a))

    drawHeu12(a,sp, x,y)
    if fsol_name != "":
        # format fsol:
        # [<len opt>]
        # [<len mst>]
        # <num sp>
        # <x> <y>
        # ...
        # <x> <y>
        optimalsp = []
        sol_file = open(fsol_name, "r")
        if fsol_has_lenghts:
            len_opt = float(sol_file.readline().strip())
            len_mst = float(sol_file.readline().strip())
        num_sp = int(sol_file.readline().strip())
        for i in range(num_sp):
            a1, a2 = sol_file.readline().strip().split(" ")
            optimalsp.append((float(a1),float(a2)))
        sol_file.close()
        plt.clf()
        drawOpt(a,sp,optimalsp)


def printMatrix(x):
    for i in range(len(x)):
        for j in range(len(x[i])):
            print('%.0f' % x[i][j], end=" ")
        print("")
    print("\n\n")

def printMatrixSideBySide(x, y) :
    for i in range(len(x)):
        for j in range(len(x[i])):
            print('%.0f' % x[i][j], end=" ")
        print(" |", end = " ")
        for j in range(len(y[i])):
            print('%.0f' % y[i][j], end=" ")
        print()
def printPoints(x):
    for p in x:
        print('(%-.7f, %-.7f)' % (p[0],p[1]), sep="\n")

def checkIsTree(x, a, sp):
    lenght = 0
    for t in range(len(a)):
        found = False
        for s in range(len(sp)):
            if x[s][t] > 0.9:
                if found:
                    print("Terminal %d linked to more sps\n" % t)
                    printMatrix(x)
                    return 0
                lenght += math.dist(a[t], sp[s])
                found = True
        if not found:
            print("Terminal %d not linked to any sp\n" % t)
            return 0
    return lenght
def len_emst(a):
    counter = itertools.count()
    REMOVED = 'removed'
    res = []
    info = {}
    queue = []
    for p in a:
        entry = [math.inf, next(counter), p, None]
        info[p] = entry
        queue.append(entry)

    while info:
        d, _, x, v = heapq.heappop(queue)
        if x is not REMOVED:
            del info[x]
            if not math.isinf(d):
                res.append((x,v))
            for _, y in info.items():
                dist = math.dist(x, y[2])
                if dist < y[0]:
                    changed_priority = [dist, next(counter), y[2], x]
                    info[y[2]] = changed_priority
                    y[2] = REMOVED
                    heapq.heappush(queue, changed_priority)
    lenght = 0
    for (x,y) in res:
        lenght += math.dist(x,y)
    return lenght
def len_heu_tree(x, a, sp):
    lenght = 0
    for t in range(len(a)):
        max_val = 0
        index = 0
        for s in range(len(sp)):
            if x[s][t] > max_val:
              max_val = x[s][t]
              index = s
        lenght += math.dist(a[t], sp[index])
    return lenght

def drawFrame(x,y,a,sp,numFrame):
    n = len(a)
    plt.cla()
    drawBig(sp, 'red',2)
    drawBig(a, 'black',2)
    for t in range(n):
        for s in range(n-2):
          if x[s][t] > 0:
            drawSegmentColouredBig(a[t], sp[s], 'black', 6, alpha=x[s][t])
    for s in range(n-2):
      for r in range(n-2):
        if y[s][r] >= 1:
            drawSegmentColouredBig(sp[r], sp[s], 'red', 6)

    fig = plt.gcf()
    fig.patch.set_facecolor((1,1,0.94))
    fig.set_size_inches(15,9)
    plt.axis('off')
    plt.grid(False)
    plt.savefig('/tmp/frame' + '{:03d}'.format(numFrame) + '.png', bbox_inches='tight')




EXE_NAME    =  "delaunay4.out"
FINPUT_NAME = "input.txt"  # file di input
ARGV = []
# ARGV = ["-n", "row"] # "row"
# formato file di input
#  <num punti>
#  <x> <y>
#  ...
#  <x> <y>
def is_power_of_two(n):
    while n > 1:
        if n % 2 == 1 :
            return False
        n /= 2
    return True
def write_terminals_to_file(f_name, a):
    input_file = open(f_name, "w")
    input_file.write(" %d\n" % len(a))
    for i in range(len(a)):
        input_file.write(" %.7f %.7f\n" % (a[i][0], a[i][1]))
    input_file.close()
def get_terminals(finput_name):
    a = []
    input_file = open(finput_name, "r")
    n = int(input_file.readline().strip())
    for i in range(n):
        p1, p2 = input_file.readline().strip().split(" ")
        a.append((float(p1),float(p2)))
    input_file.close()
    return a
def get_steiner_points(f_name):
    # formato file
    # <x>,<y>
    # ...
    # <x>,<y>
    sp = []
    file = open(f_name, "r")
    line = file.readline()
    if "INFINITE LOOP" in line:
        return None
    while line:
        p1, p2 = line.split(",")
        sp.append((float(p1),float(p2)))
        line = file.readline()
    file.close()
    return sp


def get_steiner_points_from_output_animation(f_name):
    sp = []
    file = open(f_name, "r")
    line = file.readline()
    begin = False
    while line:
        if begin:
            p1, p2 = line.split(",")
            sp.append((float(p1),float(p2)))
        elif line == "Steiner points:":
            begin = True        
        line = file.readline()
    file.close()
    return sp


def get_steiner_points_from_geo(f_name):
    # formato file
    # <x>,<y>
    # ...
    # <x>,<y>
    sp = []
    file = open(f_name, "r")
    line = file.readline()
    while line:
        if "@C" in line:
            _, p1, p2 =  line.strip().split("\t")
            sp.append((float(p1),float(p2)))
        line = file.readline()

    file.close()
    return sp


def print_results(lenghts, times, len_opt, len_mst):
    best_len = lenghts[0]
    worst_len = lenghts[0]
    for l in lenghts:
        if l < best_len:
            best_len = l
        elif l > worst_len:
            worst_len = l

    # stampe x determinazione NUM_ITERATIONS
    print("({:.2f}, {:.2f}, {:.2f})"
    .format((best_len/len_opt -1)*100, (statistics.mean(lenghts)/len_opt -1)*100,
    (worst_len/len_opt -1)*100), end=" &")

def print_results_pool(lenghts, times, len_opt, len_mst, print_to_file=""):
# restituisce info che servono x calcolo riga tabella riassuntiva
    # print("print_results_pool: file name =" + print_to_file)
    reductions = [(1 - l/len_mst)*100 for l in lenghts]
    best_red  = max(reductions)
    worst_red = min(reductions)
    avg_red   = statistics.mean(reductions)
    avg_time  = statistics.mean(times)
    opt_red = (1 - len_opt/len_mst)*100
    sse = 0
    for r in reductions:
        sse += math.pow(r - avg_red,2)

    if print_to_file == "":
        print("{:.2f} & {:.2f} & {:.2f} & {:.3f} &  {:.2f} & {:.2f}"
        .format(best_red, avg_red, worst_red, statistics.stdev(reductions), opt_red, avg_time, end=" &"))
    else:
        f = open(print_to_file, "a")
        print("{:.2f} & {:.2f} & {:.2f} & {:.3f} &  {:.2f} & {:.2f}"
        .format(best_red, avg_red, worst_red, statistics.stdev(reductions), opt_red, avg_time, end=" &"), file=f)
        f.close()

    return (best_red, avg_red, worst_red, sse, opt_red, avg_time)


def analyze_instance(finput_name, exe_name, num_ripetizioni):
    # nome file di input, nome programma da eseguire, per quante volte
    # restituisce ([len trovate], [tempi impiegati])
    lenghts = []
    times = []

    # get terminals
    a = get_terminals(finput_name)

    # get sp
    for k in range(num_ripetizioni):
        success = False
        while not success:
            if len(a) < 40: #40 per delaunay, 8 per lp
                time.sleep(1) # else non cambia il seed
            # exec program
            input_file = open(finput_name, "r")
            foutput = open("Output.txt", "w")
            total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
            subprocess.run(["./" + exe_name] + ARGV, stdin=input_file, stdout=foutput)
            new_total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
            time_spent = new_total_time - total_time
            input_file.close()
            foutput.close()
            sp = get_steiner_points("Output.txt")
            if not sp is None:
                success = True
            else:
                print("failed one", file=sys.stderr)
        len_sol = len_emst(a+sp)

        lenghts.append(len_sol)
        times.append(time_spent)

    return (lenghts, times)

def draw_instance(finput_name, exe_name, fsol_name="", fsol_has_lenghts=True, radius=0):
    a = get_terminals(finput_name)

    input_file  = open(finput_name, "r")
    output_file = open("Output.txt", "w")
    subprocess.run(["./" + exe_name], stdin=input_file, stdout=output_file)
    input_file.close()
    output_file.close()
    sp = get_steiner_points("Output.txt")

    plt.axis('equal')
    if radius == 0:
        plt.grid(True)

    
    drawHeu(a,sp, radius)
    if fsol_name != "":
        # format fsol:
        # [<len opt>]
        # [<len mst>]
        # <num sp>
        # <x> <y>
        # ...
        # <x> <y>
        optimalsp = []
        sol_file = open(fsol_name, "r")
        if fsol_has_lenghts:
            len_opt = float(sol_file.readline().strip())
            len_mst = float(sol_file.readline().strip())
        num_sp = int(sol_file.readline().strip())
        for i in range(num_sp):
            a1, a2 = sol_file.readline().strip().split(" ")
            optimalsp.append((float(a1),float(a2)))
        sol_file.close()
        plt.clf()
        drawOpt(a,sp,optimalsp)
        actual_red = (1- len_emst(a+sp)/len_mst)*100
        opt_red = (1 - len_opt/len_mst) * 100
        print("{:.2f} -- {:.2f}".format(actual_red, opt_red))

def draw_only_opt(finput_name, fsol_name, fsol_has_lenghts=True):
    a = get_terminals(finput_name)
    plt.axis('equal')
    plt.grid(True)

    # format fsol:
    # [<len opt>]
    # [<len mst>]
    # <num sp>
    # <x> <y>
    # ...
    # <x> <y>
    optimalsp = []
    sol_file = open(fsol_name, "r")
    if fsol_has_lenghts:
        len_opt = float(sol_file.readline().strip())
        len_mst = float(sol_file.readline().strip())
    num_sp = int(sol_file.readline().strip())
    for i in range(num_sp):
        a1, a2 = sol_file.readline().strip().split(" ")
        optimalsp.append((float(a1),float(a2)))
    sol_file.close()
    plt.clf()
    
    drawAnnotate(a, 'black')
    drawAnnotate(optimalsp, 'red')
    points = a + optimalsp
    matrix = emst(points)
    for t in range(len(points)):
      for s in range(len(points)):
        if matrix[s][t] == 1:
          drawSegment(points[t], points[s])

    fig = plt.gcf()
    fig.set_size_inches(16,6)
    plt.savefig('/tmp/opt.png', bbox_inches='tight')




def analyze_OR_library(num_terminals):
    NUM_ITERATIONS = 2 + math.ceil(math.log2(num_terminals))
    num_points = str(num_terminals)
    data_file = open("estein" + num_points + ".txt", "r")
    sol_file  = open("estein" + num_points + "opt.txt", "r")
    num_istanze = int(data_file.readline().strip())
    sol_file.readline()

    migliori = []
    medie = []
    peggiori = []
    somme_quadrati = []
    ottimi = []
    tempi = []

    for k in range(num_istanze):
        stdin_file = open(FINPUT_NAME, "w")
        # write input
        n = int(data_file.readline().strip())
        stdin_file.write(" %d\n" % n)
        for i in range(n):
            a1, a2 = data_file.readline().strip().split(" ")
            point = (float(a1),float(a2))
            stdin_file.write(" %.7f %.7f\n" % (point[0], point[1]))
        stdin_file.close()

        # get solution
        len_opt = float(sol_file.readline().strip())
        len_mst = float(sol_file.readline().strip())
        num_sp = int(sol_file.readline().strip())
        for i in range(num_sp):
            a1, a2 = sol_file.readline().strip().split(" ")

        lenghts, times = analyze_instance(FINPUT_NAME, EXE_NAME, NUM_ITERATIONS)
        print("%2d done" % k, file=sys.stderr)

        info = print_results_pool(lenghts, times, len_opt, len_mst)
        migliori.append(info[0])
        medie.append(info[1])
        peggiori.append(info[2])
        somme_quadrati.append(info[3])
        ottimi.append(info[4])
        tempi.append(info[5])

    sol_file.close()
    data_file.close()

    sse = 0
    for s in somme_quadrati:
        sse += s
    sse = math.sqrt(sse/(num_istanze*(NUM_ITERATIONS-1)))

    print("Summary row: best,avg,worst, stdev, ottimo, time")
    print("({:.2f} & {:.2f} & {:.2f} & {:.3f} & {:.2f} & {:.2f})"
    .format(statistics.mean(migliori), statistics.mean(medie),statistics.mean(peggiori),
    sse, statistics.mean(ottimi), statistics.mean(tempi)))

def generate_grid(r, c, len_side):
    a = []
    left_corner_X = 0
    left_corner_Y = 0
    currX = left_corner_X
    for i in range(c):
        currY = left_corner_Y
        for j in range(r):
            a.append((currX, currY))
            currY += len_side
        currX += len_side
    return a
def DEPRECATED_analyze_ladder(exe_name, num_ripetizioni, m):
    # nome programma da eseguire, per quante volte
    # m = len ladder (griglia 2xm)
    # restituisce ([len trovate], [tempi impiegati])
    lenghts = []
    times = []

    # set terminals
    a = []
    height = 1  #
    base = 0    # y delle linee lungo cui disposti i punti
    currX = 0
    for i in range(m):
        a.append((currX, base))
        a.append((currX, base+height))
        currX += height
    # write input
    n = 2*m
    write_terminals_to_file(FINPUT_NAME, a)
    
    lenghts, times = analyze_instance(FINPUT_NAME, EXE_NAME, num_ripetizioni)
    len_mst = 2*m - 1
    if m%2 == 0:
        len_opt = m * (1 + math.sqrt(3)/2) - 1
    else:
        len_opt = math.sqrt(math.pow(m * (1 + math.sqrt(3)/2) - 1, 2) + 0.25)
    print_results(lenghts, times, len_opt, len_mst)

def analyze_grid(exe_name, r,c):
    num_ripetizioni = 2 + math.ceil(math.log2(r*c))
    # nome programma da eseguire, per quante volte
    # griglia rxc
    # restituisce ([len trovate], [tempi impiegati])
    lenghts = []
    times = []

    # set terminals
    a = generate_grid(r,c, 1)
    # write input
    n = r*c
    write_terminals_to_file(FINPUT_NAME, a)    

    lenghts, times = analyze_instance(FINPUT_NAME, exe_name, num_ripetizioni)

    len_mst = r*c - 1 # se len_side = 1

    # len opt
    if r == 2: # ladder
        if c%2 == 0:
            len_opt = c * (1 + math.sqrt(3)/2) - 1
        else:
            len_opt = math.sqrt(math.pow(c * (1 + math.sqrt(3)/2) - 1, 2) + 0.25)

    elif r == c: # griglie quadrate
        rho = (1 + math.sqrt(3)) / 3
        len_X = 1 + math.sqrt(3)
        len_Y = math.sqrt(2 + math.sqrt(3))
        len_A2 = math.sqrt(11 + 6*math.sqrt(3))
        len_A4 = math.sqrt(35 + 20*math.sqrt(3))

        if r == 6:
            len_opt = (math.pow(r,2)-3)*rho + len_Y
        elif is_power_of_two(r):
            k = math.log2(r)
            len_opt = (math.pow(4,k)-1) / 3 * len_X
        else:
            match r % 6:
                case 0:
                    len_opt = (math.pow(r,2) - 4)*rho + len_A2
                case 1 | 5:
                    len_opt = (math.pow(r,2) - 4)*rho + 3
                case 2 | 4:
                    len_opt = (math.pow(r,2) - 10)*rho + len_A4
                case 3:
                    len_opt = (math.pow(r,2) - 3)*rho + 2

    print_results_pool(lenghts, times, len_opt, len_mst)

def wrapper_cocircular(exe_name, num_istanze, n, r, sigma, print_to_file=""):
    NUM_ITERATIONS = 2 + math.ceil(math.log2(n))
    migliori = []
    medie = []
    peggiori = []
    somme_quadrati = []
    ottimi = []
    tempi = []

    for k in range(num_istanze):
        if print_to_file == "":
            info = analyze_cocircular(exe_name, NUM_ITERATIONS, n, r, sigma)
        else:
            info = analyze_cocircular(exe_name, NUM_ITERATIONS, n, r, sigma, print_to_file)
        migliori.append(info[0])
        medie.append(info[1])
        peggiori.append(info[2])
        somme_quadrati.append(info[3])
        ottimi.append(info[4])
        tempi.append(info[5])
    
    sse = 0
    for s in somme_quadrati:
        sse += s
    sse = math.sqrt(sse/(num_istanze*(NUM_ITERATIONS-1)))

    if print_to_file == "":
        print("Summary row: best,avg,worst, stdev, ottimo, time")
        print("({:.2f} & {:.2f} & {:.2f} & {:.3f} & {:.2f} & {:.2f})"
        .format(statistics.mean(migliori), statistics.mean(medie),statistics.mean(peggiori),
        sse, statistics.mean(ottimi), statistics.mean(tempi)))
    else:
        f = open(print_to_file, "a")
        print("Summary row: best,avg,worst, stdev, ottimo, time", file=f)
        print("({:.2f} & {:.2f} & {:.2f} & {:.3f} & {:.2f} & {:.2f})"
        .format(statistics.mean(migliori), statistics.mean(medie),statistics.mean(peggiori),
        sse, statistics.mean(ottimi), statistics.mean(tempi)), file=f)
        f.close()
    

def analyze_cocircular(exe_name, num_ripetizioni, n, r, sigma, print_to_file=""):
    # nome programma da eseguire, per quante volte, num terminali
    # restituisce ([len trovate], [tempi impiegati])
    #
    # terminali generati in coordinate polari (ρ,θ):
    # - θ uniforme tra 0 e 2pi
    # - ρ distribuito come gaussiana di media r e stdev sigma
    lenghts = []
    times = []

    # set terminals
    random.seed()
    a = []
    for i in range(n):
        theta = random.uniform(0, 2*math.pi)
        rho = random.gauss(r, sigma)
        a.append((rho*math.cos(theta), rho*math.sin(theta)))
    # write input
    write_terminals_to_file(FINPUT_NAME, a)
    
    lenghts, times = analyze_instance(FINPUT_NAME, exe_name, num_ripetizioni)
    
    len_mst = len_emst(a)

    # exec geosteiner to get optimum
    PATH_GEO = "/home/ror/Downloads/geosteiner-5.3/"
    input_file = open(FINPUT_NAME, "r")
    input_geosteiner = open("input_geosteiner.txt", "w")
    input_file.readline() # skip number of terminals
    line = input_file.readline()
    while line:
        input_geosteiner.write(line)
        line = input_file.readline()
    input_file.close()
    input_geosteiner.close()

    foutput = open("Output.txt", "w")
    total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    p1 = subprocess.run(["cat", "input_geosteiner.txt"], stdout=subprocess.PIPE)
    p2 = subprocess.run([PATH_GEO + "efst"], input=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.run([PATH_GEO + "bb"], input=p2.stdout, stdout=foutput)
    new_total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    time_spent = new_total_time - total_time
    foutput.close()

    sp = get_steiner_points_from_geo("Output.txt")

    len_opt = len_emst(a+sp)
    if print_to_file == "":
        return print_results_pool(lenghts, times, len_opt, len_mst)
    else:
        return print_results_pool(lenghts, times, len_opt, len_mst, print_to_file)

def analyze_big(exe_name, num_ripetizioni, n):
    lenghts = []
    times = []

    # set terminals
    random.seed()
    a = []
    for i in range(n):
        x = random.uniform(0,1)
        y = random.uniform(0,1)
        a.append((x, y))
    # write input
    write_terminals_to_file(FINPUT_NAME, a)
    
    lenghts, times = analyze_instance(FINPUT_NAME, exe_name, num_ripetizioni)
    
    len_mst = len_emst(a)

    # exec geosteiner to get optimum
    PATH_GEO = "/home/ror/Downloads/geosteiner-5.3/"
    input_file = open(FINPUT_NAME, "r")
    input_geosteiner = open("input_geosteiner.txt", "w")
    input_file.readline() # skip number of terminals
    line = input_file.readline()
    while line:
        input_geosteiner.write(line)
        line = input_file.readline()
    input_file.close()
    input_geosteiner.close()

    foutput = open("Output.txt", "w")
    total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    p1 = subprocess.run(["cat", "input_geosteiner.txt"], stdout=subprocess.PIPE)
    p2 = subprocess.run([PATH_GEO + "efst"], input=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.run([PATH_GEO + "bb"], input=p2.stdout, stdout=foutput)
    new_total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    time_spent = new_total_time - total_time
    foutput.close()

    sp = get_steiner_points_from_geo("Output.txt")
    #sp = get_steiner_points("Output.txt")

    len_opt = len_emst(a+sp)
    if print_to_file == "":
        return print_results_pool(lenghts, times, len_opt, len_mst)

def animate(exec=False):
    if exec:
        singleRunxAnimate(FINPUT_NAME, EXE_NAME)

    a = get_terminals(FINPUT_NAME)
    n = len(a)

    x  = []
    y  = []
    sp = []
    numFrame = 0

    file = open("/tmp/Output.txt", "r")
    readingMatrix  = True

    line = file.readline()
    while line != "Steiner points:":
        if line == "\n":
            readingMatrix = False
        elif readingMatrix:
            if len(line.split("|")) == 1:
                print(line)
            xs, ys = line.split("|")
            x.append([float(c) for c in xs.strip().split(" ")])
            y.append([float(c) for c in ys.strip().split(" ")])
        else:
            for index in range(n-2):
                a1, a2 = line.split(",")
                sp.append((float(a1),float(a2)))
                line = file.readline()
            drawFrame(x,y,a,sp,numFrame)
            numFrame += 1
            x  = []
            y  = []
            sp = []
            readingMatrix = True
            # break
        line = file.readline()
    file.close()

def singleRunxAnimate(finput_name, exe_name):
    """
    exec and print len of sol find. On /tmp/Output.txt infos for creating animation.
    """
    # get terminals
    a = get_terminals(finput_name)

    # get sp
    success = False
    while not success:
        # exec program
        input_file = open(finput_name, "r")
        foutput = open("/tmp/Output.txt", "w")
        total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
        subprocess.run(["./" + exe_name] + ARGV, stdin=input_file, stdout=foutput)
        new_total_time = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
        time_spent = new_total_time - total_time
        input_file.close()
        foutput.close()
        sp = get_steiner_points_from_output_animation("/tmp/Output.txt")
        if not sp is None:
            success = True
        else:
            print("failed one", file=sys.stderr)
    len_sol = len_emst(a+sp)

    print(len_sol)


animate()