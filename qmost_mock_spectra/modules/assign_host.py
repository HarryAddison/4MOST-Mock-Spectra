import random


def _assign_random_host():
    # cumulative probability bins:
    # E, E/S0, S0, S0/a, Sa, Sab, Sb, Sbc, Sc, Scd, Sd+
    Ia = [0.08356164383561644, 0.1315068493150685, 0.20136986301369864, 0.3013698630136986, 0.34383561643835614, 0.41232876712328764, 0.5794520547945206, 0.7506849315068493, 0.9068493150684932, 0.9452054794520548, 1.0]

    rand_num = random.random()

    if rand_num < Ia[0]:
        host = "e"
    elif Ia[0] <= rand_num < Ia[1]:
        host = random.choice(["e", "s0"])
    elif Ia[1] <= rand_num < Ia[2]:
        host = "s0"
    elif Ia[2] <= rand_num < Ia[3]:
        host = random.choice(["s0", "sa"])    
    elif Ia[3] <= rand_num < Ia[4]:
        host = "sa"
    elif Ia[4] <= rand_num < Ia[5]:
        host = random.choice(["sa","sb"])
    elif Ia[5] <= rand_num < Ia[6]:
        host = "sb"
    elif Ia[6] <= rand_num < Ia[7]:
        host = random.choice(["sb", "sc"])
    elif Ia[7] <= rand_num < Ia[8]:
        host = "sc"
    elif Ia[8] <= rand_num < Ia[9]:
        host = "sc"
    elif rand_num > Ia[9]:
        host = random.choice(["sb","sc"])
    else:
        print(rand_num, ' is not recognised for some reason')

    return host


def assign_host(host_type="random"):

    if host_type == "random":
        host_type = _assign_random_host()
    return host_type