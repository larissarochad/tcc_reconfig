#--------------------------------
# Reconfiguracao
#--------------------------------

#--------------------------------
# definir os conjuntos
#--------------------------------

set Ob;
set Ol within Ob cross Ob;

#--------------------------------
# definir os parametros
#--------------------------------

# dados das barras
param Tb{Ob};     # tipo de barra (1: barra SE, 2: barra de carga)
param PD{Ob};     # potencia ativa demanda
param QD{Ob};     # potencia reativa demanda
param QC{Ob};     # capacitor existente na barra

# dados das linhas
param R{Ol};      # resistencia das linhas
param X{Ol};      # reatancia indutiva das linhas
param Imax{Ol};   # corrente maxima das linhas
param Z2{Ol};     # Z^2 = R^2 + X^2
param linha{Ol};
param l{Ol};
param cond{Ol};

param Vnom;
param Vmin;
param Vmax;

param Ke = 168;

#--------------------------------
# para a linearização Vsqr * Isqr
#--------------------------------

param S = 0;
param deltaV = (Vmax^2 - Vmin^2) / (S + 1);
var x{Ob, s in 1..S}, binary;
var Pc{Ob, s in 1..S};

#--------------------------------
# para a linearização de P^2 e Q^2
#--------------------------------

param Y = 50;
param deltaS{Ol};
param ms{Ol, y in 1..Y};
var deltaP{Ol, y in 1..Y} >= 0;
var deltaQ{Ol, y in 1..Y} >= 0;

var Pmais{Ol} >= 0;
var Pmenos{Ol} >= 0;
var Qmais{Ol} >= 0;
var Qmenos{Ol} >= 0;

#--------------------------------
# definir as variaveis
#--------------------------------

var PS{Ob};    # potencia ativa da subestacao
var QS{Ob};    # potencia reativa da subestacao

var P{Ol};     # potencia ativa na linha
var Q{Ol};     # potencia reativa na linha

var Isqr{Ol} >= 0;  # corrente nas linhas
var Vsqr{Ob};  # tensoes nas barras

var b{Ol};
var ymais{Ol}, binary;
var ymenos{Ol}, binary;

#--------------------------------
# funcao objetivo

minimize FuncaoObjetivo: (Ke * sum{(i,j) in Ol}(R[i,j] * Isqr[i,j]));

#--------------------------------
# balanco de potencia ativa

subject to BalacoPotenciaAtiva{i in Ob}:
	sum{(k,i) in Ol}(P[k,i]) - sum{(i,j) in Ol}(P[i,j] + R[i,j] * Isqr[i,j]) + PS[i] = PD[i];
	
#--------------------------------
# balanco de potencia reativa

subject to BalacoPotenciaReativa{i in Ob}:
	sum{(k,i) in Ol}(Q[k,i]) - sum{(i,j) in Ol}(Q[i,j] + X[i,j] * Isqr[i,j]) + QS[i] = QD[i];

#--------------------------------
# queda de tensao no circuito

subject to QuedaTensao{(i,j) in Ol}:
	Vsqr[i] - 2*(R[i,j] * P[i,j] + X[i,j] * Q[i,j]) - Z2[i,j] * Isqr[i,j] - Vsqr[j] - b[i,j] = 0;

#--------------------------------
# potencia aparente (kVA)

subject to PotenciaAparente{(i,j) in Ol}:
	(Vmin^2 + (1/2) * deltaV) * Isqr[i,j] + (sum{s in 1..S} (Pc[j,s])) = (sum{y in 1..Y} (ms[i,j,y] * deltaP[i,j,y])) + (sum{y in 1..Y} (ms[i,j,y] * deltaQ[i,j,y]));
	
#--------------------------------
# potencia ativa

subject to LinearAtiva1{(i,j) in Ol}:
	Pmais[i,j] - Pmenos[i,j] = P[i,j];
	
subject to LinearAtiva2{(i,j) in Ol}:
	Pmais[i,j] + Pmenos[i,j] = sum{y in 1..Y} (deltaP[i,j,y]);

subject to LinearAtiva3{(i,j) in Ol, y in 1..Y}:
	deltaP[i,j,y] <= deltaS[i,j];

#--------------------------------
# potencia reativa

subject to LinearReativa1{(i,j) in Ol}:
	Qmais[i,j] - Qmenos[i,j] = Q[i,j];

subject to LinearReativa2{(i,j) in Ol}:
	Qmais[i,j] + Qmenos[i,j] = sum{y in 1..Y} (deltaQ[i,j,y]);

subject to LinearReativa3{(i,j) in Ol, y in 1..Y}:
	deltaQ[i,j,y] <= deltaS[i,j];
	
#--------------------------------
# limite das tensoes

subject to LimiteTensao{i in Ob}:
	Vmin^2 <= Vsqr[i] <= Vmax^2;

#--------------------------------
# limite de corrente

subject to LimiteCorrente{(i,j) in Ol}:
	Isqr[i,j] <= Imax[i,j]^2 * (ymais[i,j] + ymenos[i,j]);

#--------------------------------
# limite das tensoes

subject to LimiteTensao2_1{j in Ob}:
	Vmin^2 + (sum{s in 1..S} (x[j,s] * deltaV)) <= Vsqr[j];
	
subject to LimiteTensao2_2{j in Ob}:
	Vsqr[j] <= Vmin^2 + deltaV + (sum{s in 1..S} (x[j,s] * deltaV));

#--------------------------------
# valores das correcoes de potencia

subject to PotenciaCorrecao1_1{(i,j) in Ol, s in 1..S}:
	0 <= deltaV * Isqr[i,j] - Pc[j,s];

subject to PotenciaCorrecao1_2{(i,j) in Ol, s in 1..S}:
	deltaV * Isqr[i,j] - Pc[j,s] <= deltaV * Imax[i,j]^2 * (1 - x[j,s]);

#--------------------------------
# valores das correcoes de potencia

subject to PotenciaCorrecao2_1{(i,j) in Ol, s in 1..S}:
	0 <= Pc[j,s];

subject to PotenciaCorrecao2_2{(i,j) in Ol, s in 1..S}:
	Pc[j,s] <= deltaV * Imax[i,j]^2 * x[j,s];

#--------------------------------
# variavel binaria

subject to VariavelBinaria{j in Ob, s in 2..S}:
	x[j,s] <= x[j,s-1];

#--------------------------------
# 

subject to Equation34{(i,j) in Ol}:
	Pmais[i,j] <= Vmax * Imax[i,j] * ymais[i,j];

#--------------------------------
# 

subject to Equation35{(i,j) in Ol}:
	Pmenos[i,j] <= Vmax * Imax[i,j] * ymenos[i,j];

#--------------------------------
# 

subject to Equation36_1{(i,j) in Ol}:
	-(Vmax * Imax[i,j] * (ymais[i,j] + ymenos[i,j])) <= Q[i,j];

subject to Equation36_2{(i,j) in Ol}:
	Q[i,j] <= Vmax * Imax[i,j] * (ymais[i,j] + ymenos[i,j]);

#--------------------------------
# 

subject to Equation37_1{(i,j) in Ol}:
	-((Vmax^2 - Vmin^2) * (1 - (ymais[i,j] + ymenos[i,j]))) <= b[i,j];

subject to Equation37_2{(i,j) in Ol}:
	b[i,j] <= (Vmax^2 - Vmin^2) * (1 - (ymais[i,j] + ymenos[i,j]));

#--------------------------------
# 

subject to Equation38:
	sum{(i,j) in Ol}(ymais[i,j] + ymenos[i,j]) = card(Ob) - 1;

#--------------------------------
# 

subject to Equation39{(i,j) in Ol}:
	(ymais[i,j] + ymenos[i,j]) <= 1;

#--------------------------------
# 

subject to jabra1{i in Ob: Tb[i] == 1}:
    sum{(i,j) in Ol}(ymais[i,j]) + sum{(j,i) in Ol}(ymenos[j,i]) >= 1;

#--------------------------------
# 
subject to jabra2{i in Ob: Tb[i] == 0}:
    sum{(i,j) in Ol}(ymenos[i,j]) + sum{(j,i) in Ol}(ymais[j,i]) = 1;