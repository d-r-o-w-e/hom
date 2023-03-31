import math
from itertools import product

# ~~~ CLASSES ~~~

class AbGp:
    def __init__(self, summands=[]):
        if not summands:
            self.summands = ["0"]
        else:
            self.summands = summands
            
    def __str__(self):
        return "+".join(list(sorted(self.summands)))
    
    def __repr__(self):
        return "+".join(list(sorted(self.summands)))
    
class Space:
    def __init__(self, root, children=[]):
        self.root = root
        self.children = children
        
    def __str__(self):
        if len(self.children) == 0:
            return self.root
        elif len(self.children) == 1:
            return self.root + "(" + str(self.children[0]) + ")"
        elif len(self.children) == 2:
            return str(self.children[0]) + self.root + str(self.children[1])
        else:
            raise IndexError("Too many child nodes.")
    
    def __repr__(self):
        if len(self.children) == 0:
            return "Space(" + self.root + ")"
        else:
            return "Space(" + self.root + ", " + repr(self.children) + ")"
    
# ~~~ AB GP OPERATIONS ~~~

def dsum(ab1, ab2):
    return AbGp(list(filter(lambda x: x!="0", ab1.summands + ab2.summands)))

def tens(ab1, ab2):
        # define tensor of single groups
        def tens1(x, y):
            if x == "R" or y == "R": # hack to get around implementing R, since i've already implemented Q
                x_new = "Q" if x == "R" else x
                y_new = "Q" if y == "R" else y
                return "R" if tens1(x_new, y_new) == "Q" else tens1(x_new, y_new)
            elif x == "Z": # if either is Z, return the other
                if y != "0":
                    return y
                elif y == "0":
                    return "0"
            elif x[0] == "Z":
                if y == "Z":
                    return x
                elif y[0] == "Z":
                    if math.gcd(int(x[1:]), int(y[1:])) != 1:
                        return "Z" + str(math.gcd(int(x[1:]), int(y[1:])))
                    else:
                        return "0"
                elif y == "Q":
                    return "0"
                elif y == "0":
                    return "0"
                else:
                    raise NotImplementedError("Unknown group encountered")
            elif x == "Q":
                if y == "Z":
                    return "Q"
                elif y[0] == "Z":
                    return "0"
                elif y == "Q":
                    return "Q"
                elif y == "0":
                    return "0"
            elif x == "0":
                return "0"
            else:
                raise NotImplementedError("Unknown group encountered")
        
        # next, find the direct sum of tensoring each summand in each part, removing 0s
        tens_prod = list(filter(lambda x: x!="0", [tens1(ab1.summands[i], ab2.summands[j]) for i, j in product(range(len(ab1.summands)), range(len(ab2.summands)))]))
        
        # return the resulting abgp
        return AbGp(tens_prod)
    
def Tor(ab1, ab2):
    # Tor functor definitions
    def Tor1(x, y):
        if x == "R" or y == "R": # hack to get around implementing R, since i've already implemented Q
                x_new = "Q" if x == "R" else x
                y_new = "Q" if y == "R" else y
                return "R" if Tor1(x_new, y_new) == "Q" else Tor1(x_new, y_new)
        elif x == "Z" or x == "0" or y == "0":
            return "0"
        elif x[0] == "Z":
            if y == "Z" or y == "0" or y == "Q":
                return "0"
            elif y[0] == "Z":
                if math.gcd(int(x[1:]), int(y[1:])) != 1:
                    return "Z" + str(math.gcd(int(x[1:]), int(y[1:])))
                else:
                    return "0"
            else:
                raise NotImplementedError("Unknown group encountered")
        elif x == "Q":
            if y == "Z" or y == "0" or y[0] == "0":
                return "0"
            elif y == "Q":
                raise NotImplementedError("Not sure what Tor(Q, Q) is")
        elif x == "R":
            if y == "Z" or y == "0" or y[0] == "0":
                return "0"
            elif y == "Q":
                raise NotImplementedError("Not sure what Tor(Q, Q) is")
        else:
            raise NotImplementedError("Unknown group encountered")
    
    Tor_prod = list(filter(lambda x: x!="0", [Tor1(ab1.summands[i], ab2.summands[j]) for i, j in product(range(len(ab1.summands)), range(len(ab2.summands)))]))
    
    return AbGp(Tor_prod)

def Ext(ab1, ab2):
    def Ext1(x, y):
        if x == "R" or y == "R": # hack to get around implementing R, since i've already implemented Q
                x_new = "Q" if x == "R" else x
                y_new = "Q" if y == "R" else y
                return "R" if Ext1(x_new, y_new) == "Q" else Ext1(x_new, y_new)
        elif x == "Z" or x == "0" or y == "0":
            return "0"
        elif x[0] == "Z":
            if y == "Z":
                return x
            elif y[0] == "Z":
                if math.gcd(int(x[1:]), int(y[1:])) != 1:
                    return "Z" + str(math.gcd(int(x[1:]), int(y[1:])))
                else:
                    return "0"
            elif y == "Q":
                return "0"
        elif y == "Q":
            return "0"
        else:
            raise NotImplementedError("Unknown group encountered")
    
    # doesn't know how to handle second ext part.  thats fine for our purposes
    Ext_prod = list(filter(lambda x: x!="0", [Ext1(ab1.summands[i], ab2.summands[0]) for i in range(len(ab1.summands))]))
    
    return AbGp(Ext_prod)

def Hom(ab1, ab2):
    def Hom1(x, y):
        if x == "R" or y == "R": # hack to get around implementing R, since i've already implemented Q
                x_new = "Q" if x == "R" else x
                y_new = "Q" if y == "R" else y
                return "R" if Hom1(x_new, y_new) == "Q" else Hom1(x_new, y_new)
        elif x == "Z":
            return y
        elif x == "0" or y == "0":
            return "0"
        elif x[0] == "Z":
            if y == "Z":
                return "0"
            elif y[0] == "Z":
                if math.gcd(int(x[1:]), int(y[1:])) != 1:
                    return "Z" + str(math.gcd(int(x[1:]), int(y[1:])))
                else:
                    return "0"
            elif y == "Q":
                return "0"
            else:
                raise NotImplementedError("Unknown group encountered")
        elif x == "Q":
            if y == "Q":
                return "Q"
            else:
                raise NotImplementedError("Unknown group encountered")
        else:
                raise NotImplementedError("Unknown group encountered")
    
    Hom_prod = list(filter(lambda x: x!="0", [Hom1(ab1.summands[i], ab2.summands[j]) for i, j in product(range(len(ab1.summands)), range(len(ab2.summands)))]))
    
    return AbGp(Hom_prod)

# ~~~ HOMOLOGY COMPUTATION ~~~

def remove_trailing_zeros(hom):
    # take a homology and remove trailing 0 groups
    while len(hom[-1].summands) == 1 and hom[-1].summands[0] == "0":
        hom = hom[:-1]
    return hom

def homology_w_coeffs(int_hom, coeff):
    # apply UCT to homology groups  
    # add a zero on the end so there are no issues, then remove it after
    int_hom += [AbGp(["0"])]
    return remove_trailing_zeros([dsum(tens(int_hom[i], coeff), Tor(int_hom[i-1], coeff)) for i in range(len(int_hom))])

def kunneth(hom1, hom2):
    # apply the kunneth theorem to 2 homology groups
    
    # pad the first and second groups to the same length:
    while len(hom1) != len(hom2):
        if len(hom1) > len(hom2):
            hom2 += [AbGp(["0"])]
        else:
            hom1 += [AbGp(["0"])]
    
    # main kunneth loop
    N = len(hom1)
    
    # add an extra 0 to make it work
    hom1 += [AbGp(["0"])]
    hom2 += [AbGp(["0"])]
    
    curr_hom = []
    for n in range(N+1): # go to N+1 to capture tor of the previous group
        # construct the abelian group coming from the tensor product of the two groups
        tens_ab = AbGp(["0"])
        p = n     # the two indices, which must always equal n
        q = n - p
        while q <= n:
            tens_ab = dsum(tens_ab, tens(hom1[p], hom2[q]))
            p -= 1 # traverse the diagonal
            q += 1
        # construct the abelian group coming from the Tor term
        tor_ab = AbGp(["0"])
        if n != 0:
            p = n-1
            q = (n-1) - p
            while q <= n-1:
                tor_ab = dsum(tor_ab, Tor(hom1[p], hom2[q]))
                p -= 1
                q += 1
        # direct sum the 2 terms, since the kunneth sequence splits
        curr_hom += [dsum(tens_ab, tor_ab)]
    
    return remove_trailing_zeros(curr_hom)

def wedge(hom1, hom2):
    # apply the wedge product to homology groups
    # this just adds everything except in dimension 0, where it gives only 1 copy of Z
    
    # pad to same length
    while len(hom1) != len(hom2):
        if len(hom1) > len(hom2):
            hom2 += [AbGp(["0"])]
        else:
            hom1 += [AbGp(["0"])]
    
    # add up like dimensions except for the first
    curr_hom = [AbGp(["Z"])] 
    if len(hom1) > 1:
        for i in range(1, len(hom1)):
            curr_hom += [dsum(hom1[i], hom2[i])]
    
    return remove_trailing_zeros(curr_hom)

def susp(hom1):
    # apply suspension to the homology groups
    # moves everything up a dimension, except keeps one copy of Z in the zeroth dimension
    # assume H0 = sum of Zs
    
    hom_susp = [AbGp(["Z"])] # k = 0
    
    # for k = 1, remove one copy of 0 from the first dimension
    if len(hom1[0].summands) == 1:
        hom_susp += [AbGp(["0"])]
    else:
        hom_susp += [AbGp(hom1[0].summands[:-1])]
    
    # for larger dimensions, just tack on the extra stuff
    if len(hom1) >= 2:
        hom_susp += hom1[1:]
    
    return remove_trailing_zeros(hom_susp)

def disj_U(hom1, hom2):
    # apply disjoint union to the homology groups
    # ... just adds up in every dimension
    
    # pad to same length
    while len(hom1) != len(hom2):
        if len(hom1) > len(hom2):
            hom2 += [AbGp(["0"])]
        else:
            hom1 += [AbGp(["0"])]
    
    # add up like dimensions except for the first
    curr_hom = [] 
    for i in range(len(hom1)):
        curr_hom += [dsum(hom1[i], hom2[i])]
    
    return remove_trailing_zeros(curr_hom)

def homology(sp1, coeff=AbGp(["Z"])):
    # given a space and an abelian group, return the homology as a list of abelian groups
    # supports spaces which have roots either a base space, or roots as "x", "v", "_", "$"
    
    # if there are strange coefficients, apply UCT for homology and return
    if coeff.summands[0] != "Z" or len(coeff.summands) != 1:
        return homology_w_coeffs(homology(sp1, AbGp(["Z"])), coeff)
    
    # base spaces occur when the space has no children
    if len(sp1.children) == 0:
        if sp1.root[0] == "S":
            # the sphere Sn has Z in dims 0 and n
            n = int(sp1.root[1:])
            curr_hom = []
            if n == 0:
                curr_hom = [AbGp(["Z", "Z"])]
            else:
                for i in range(n + 1):
                    if i == 0 or i == n:
                        curr_hom += [AbGp(["Z"])]
                    else:
                        curr_hom += [AbGp(["0"])]
            return remove_trailing_zeros(curr_hom)
        elif sp1.root[0:2] == "RP":
            # RPn has Z in dims 0 and odd n, and Z2 otherwise if odd, and 0 otherwise
            n = int(sp1.root[2:])
            curr_hom = []
            for i in range(n + 1):
                if i == 0 or (i == n and n % 2 == 1):
                    curr_hom += [AbGp(["Z"])]
                elif i % 2 == 1 and i < n:
                    curr_hom += [AbGp(["Z2"])]
                else:
                    curr_hom += [AbGp(["0"])]
            return remove_trailing_zeros(curr_hom)
        elif sp1.root[0:2] == "CP":
            # CPn has Z in odd dimensions up to 2n
            n = int(sp1.root[2:])
            curr_hom = []
            for i in range(2*n + 1):
                if i % 2 == 0:
                    curr_hom += [AbGp(["Z"])]
                else:
                    curr_hom += [AbGp(["0"])]
            return remove_trailing_zeros(curr_hom)
        elif sp1.root[0] == "D" or sp1.root == "*":
            # the disc/point in any dimension has only Z in the zeroth position
            return [AbGp(["Z"])]
        else:
            raise NotImplementedError("Encountered unknown space")
    
    # if the space does have children, it is a combination of base spaces.  use helper functions to compute those:
    if sp1.root == "x":
        return kunneth(homology(sp1.children[0], AbGp(["Z"])), homology(sp1.children[1], AbGp(["Z"])))
    elif sp1.root == "v":
        return wedge(homology(sp1.children[0], AbGp(["Z"])), homology(sp1.children[1], AbGp(["Z"])))
    elif sp1.root == "_":
        return disj_U(homology(sp1.children[0], AbGp(["Z"])), homology(sp1.children[1], AbGp(["Z"])))
    elif sp1.root == "$":
        return susp(homology(sp1.children[0], AbGp(["Z"])))
    else:
        raise NotImplementedError("Encountered unknown combinator")

def homology_str(sp, coeff=AbGp(["Z"])):
    # prettier output of homology
    hom = homology(sp, coeff)
    table_header = "  dim      |"
    def pad_to_max(s, i):
        max_len = max(len(str(i)), len(str(hom[i])))
        return s + (" "*(max_len-len(s)))
    for i in range(len(hom)):
        table_header += " " + pad_to_max(str(i), i) + " "
    hom_part = "  homology |"
    for i in range(len(hom)):
        hom_part += " " + pad_to_max(str(hom[i]), i) + " "
    return table_header + "\n" + hom_part

def cohomology(sp1, coeff=AbGp(["Z"])):
    # use UCT for cohomology to compute the cohomology
    
    # first compute the int homology
    int_hom = homology(sp1, AbGp(["Z"]))
    # add zero on the end to prevent issues; remove it after
    int_hom += [AbGp(["0"])]
    
    # then apply the UCT for cohomology:
    int_cohom = remove_trailing_zeros([dsum(Ext(int_hom[i-1], AbGp(["Z"])), Hom(int_hom[i], AbGp(["Z"]))) for i in range(len(int_hom))])
    
    # then apply tens UCT for cohomology, if necessary:
    cohom = []
    if coeff.summands[0] != "Z" or len(coeff.summands) != 1:
        int_cohom += [AbGp(["0"])]
        cohom = remove_trailing_zeros([dsum(tens(int_cohom[i], coeff), Tor(int_cohom[i+1], coeff)) for i in range(len(int_cohom) - 1)])
    else:
        cohom = int_cohom
    
    return cohom

def cohomology_str(sp, coeff=AbGp(["Z"])):
    # prettier output of homology
    hom = cohomology(sp, coeff)
    table_header = "  dim        |"
    def pad_to_max(s, i):
        max_len = max(len(str(i)), len(str(hom[i])))
        return s + (" "*(max_len-len(s)))
    for i in range(len(hom)):
        table_header += " " + pad_to_max(str(i), i) + " "
    hom_part = "  cohomology |"
    for i in range(len(hom)):
        hom_part += " " + pad_to_max(str(hom[i]), i) + " "
    return table_header + "\n" + hom_part

# ~~~ STRING PARSING ~~~

def tokenize(s):
    s = [*s.replace(" ", "")] + [0]# strip whitespace, convert to list, add end character
    tokens = []
    prev_i = 0
    running_str = ""
    for i in range(len(s)):
        running_str = "".join(s[prev_i:i])
        if s[i] in ["$", "v", "x", "_", "(", ")", 0]:
            if len(running_str) >= 1:
                tokens += [running_str, s[i]] if s[i] != 0 else [running_str] # don't include the end character
            else:
                tokens += [s[i]]
            running_str = ""
            prev_i = i + 1
    return tokens

def build_space(tokens):
    # recursive descent algorithm, following https://www.engr.mun.ca/~theo/Misc/exp_parsing.htm
    tokens += [0]  # use 0 to represent end token
    i = 0
    
    # operators
    U = {"$"}
    B = {"x", "v", "_"}
    sentinel = 1
    
    # basic subroutines
    def next():
        if i < len(tokens):
            return tokens[i]
        else:
            return tokens[i-1]
    def consume():
        nonlocal i
        if next() == 0:
            return 0
        else:
            t = tokens[i]
            i += 1
            return t
    def error():
        raise IndexError("Parsing error.")
    def expect(tok):
        if next() == tok:
            consume()
        else:
            error()   
    
    operators = []
    operands = []
    
    def EParser():
        nonlocal operators
        nonlocal operands
        operators += [1]
        E()
        expect(0)
        return operands[-1]
    
    def E():
        nonlocal operators
        nonlocal operands
        P()
        while next() in B:
            pushOp(next())
            consume()
            P()
        while operators[-1] != 1:
            popOp()
          
    def P():
        nonlocal operators
        nonlocal operands
        if not (next() in B or next() in U or next() == 0 or next() == 1 or next() == "(" or next() == ")"):
            operands += [Space(next())]
            consume()
        elif next() == "(":
            consume()
            operators += [1]
            E()
            expect(")")
            operators = operators[:-1]
        elif next() in U:
            pushOp(next())
            consume()
            P()
        else:
            error()
    
    def popOp():
        nonlocal operators
        nonlocal operands
        if operators[-1] in B:
            t1 = operands[-1]
            operands = operands[:-1]
            t0 = operands[-1]
            operands = operands[:-1]
            operands += [Space(operators[-1], [t0, t1])]
            operators = operators[:-1]
        else:
            oper = operators[-1]
            opnd = operands[-1]
            operators = operators[:-1]
            operands = operands[:-1]
            operands += [Space(oper, [opnd])]
          
    def pushOp(op):
        nonlocal operators
        nonlocal operands
        op_int = {"$": 4, "x": 3, "v": 2, "_": 1, sentinel: 0, }
        while op_int[operators[-1]] > op_int[op]:
            popOp()
        operators += [op]
        
    return EParser()

# ~~~ SOME SPACES FOR TESTING ~~~

RP2 = Space("RP2")
RP10 = Space("RP10")
S2 = Space("S2")
T2 = Space("x", [Space("S1"), Space("S1")])
S2vS1 = Space("v", [Space("S1"), Space("S2")])
S1xD2 = Space("x", [Space("S1"), Space("D2")])
SS1 = Space("$", [Space("S1")])
SSS1 = Space("$", [SS1])
S2_by_susp = Space("$", [Space("$", [Space("_", [Space("*"), Space("*")])])])
point_susp = Space("$", [Space("*")])

if __name__=="__main__":    
    while 1:
        print("Input H(*) <space description>; coeff:")
        inp = input()
        t = inp.split(";")
        if len(t) == 2:
            l, r = t
            if l[:2] != "H " and l[:3] != "H* ":
                print("Choose H or H*")
            else:
                coeffs = AbGp(r.replace(" ", "").split("+"))
                hom = True
                if l[:2] == "H ":
                    space = build_space(tokenize(l[2:]))
                    print(homology_str(space, coeffs))
                elif l[:3] == "H* ":
                    space = build_space(tokenize(l[3:]))
                    print(cohomology_str(space, coeffs))
        elif len(t) == 1:
            l = t[0]
            if l[:2] != "H " and l[:3] != "H* ":
                print("Choose H or H*")
            else:
                coeffs = AbGp(["Z"])
                hom = True
                if l[:2] == "H ":
                    space = build_space(tokenize(l[2:]))
                    print(homology_str(space, coeffs))
                elif l[:3] == "H* ":
                    space = build_space(tokenize(l[3:]))
                    print(cohomology_str(space, coeffs))
        else:
            print("Could not parse command.")
        print("")