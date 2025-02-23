include("HODLRStruct.jl")

# building a holder matrix 
# in theory we should do thisw with func constructHOLDR
H1 = HODLR(LR1= LR, LR2 = LR, H1 = array, H2 = array )
H2 = HODLR(LR1= newLR, LR2 = LR, H1 = array, H2 = array )
L1 = LR(U, V)
L2 = LR(U1,V2)

Big_H = HODLRMat(LR1 = L1, LR2 = L2, H1=H1, H2=H2)


# lets do operation 
# have a vector x

y = Big_H*x
y = Big_H@ mat
y = Big_H + mat
😶