import os
import re
import matplotlib.pyplot as plt
import numpy as np
f = open("profile1180.data","r")
lines = f.readlines()
count = 0
pattern_digit = r"\s+([+-]?\d+.?\d*[eE]?[-+]?\d*)"
M = []
C12 = []
Rho = []
M_S = 1.989e33
R_S = 696340e5
R_P = []
Result = []
Result_Dict = []
for i in range(len(lines)):
    if(i<=4):
        continue
    line = lines[i]
    pattern2 = re.compile(r"[^\s]+")
    if(i==5):
        Result_Dict = pattern2.findall(line)
        Result = [[] for i in range(len(Result_Dict)+1)]
        Result[0] = Result_Dict
        Result[0] = [0] + Result[0]
        continue
    pattern1 = re.compile(pattern_digit)
    matches = pattern1.finditer(line)
    output = ""
    if(matches):
        element_index = 0
        for match in matches:
            element_index += 1
            element = match.group(1)
            Result[element_index].append(float(element))
R_P = np.pow(10,np.array(Result[3]))*R_S
plt.xscale("log")
start = 71
end = 93
cutting = next((i for i, x in enumerate(R_P) if x < 1.5e10), None)
for i in range(start,end):
    plt.plot(R_P[cutting:],Result[i][cutting:],label=Result[0][i])
legend = plt.legend()
for i in range(end-start):
    legend_colors = [handle.get_color() for handle in legend.legend_handles]

    colors_i = legend_colors[i]
    for x in range(len(Result[start+i])-1,cutting-1,-1):
        if Result[start + i][x]>0.2*np.max(Result[start + i][cutting:]):
            plt.text(R_P[x], Result[start+i][x], Result[0][start+i],color=colors_i)
            break
plt.title("From Progenitor")
plt.show()
for i in range(27,31):
    plt.plot(R_P,Result[i],label=Result[0][i])
    plt.xscale("log")
    plt.legend()
    plt.show()
#ctrl +/ 多行注释/解除
# isFirst = False
# R = []
# count = 0
# lm = 0
# lr = 0
# M = list(reversed(M))
# Rho = list(reversed(Rho))
# dm = 0
# dr = 0
# for m in M:
#     print(m)
#     rho = Rho[count]
#     count += 1
#     if(not isFirst):
#         isFirst = True
#         r = np.pow(3*m/(4*np.pi*rho),1/3)
#         lr = r
#         lm = m
#         R.append(r)
#         print(rho,r,m)
#         print("m=", m, "lm=", lm, "dm=", dm)
#         print("r=", r / 1e5, "km ", "lr=", lr / 1e5, "km ", "dr=", dr / 1e5, "km ")
#         continue
#     dm = m-lm
#     print("m=",m,"lm=",lm,"dm=",dm)
#     dr = dm/(4*np.pi*(lr**2)*rho)
#     r = lr+dr
#     print("r=",r/1e5,"km ","lr=",lr/1e5,"km ","dr=",dr/1e5,"km ")
#     R.append(r)
#     lm = m
#     lr = r
#     #os.system("pause")
# R = list(reversed(R))
# plt.xscale("log")
# plt.plot(R,C12)
# plt.show()
