using CSV
using DataFrames
"""
function get_C, input area A and out put C from a csv file. 
csv file dataformat is A,C
"""
function get_C(A::Float64)
    #read the csv file
    filename = "AtoC.csv" #Make sure you got csv in the name
    data = CSV.read(filename,DataFrame)
    #get the column of A
    A_data = data[:,1]
    #get the column of C
    C_data = data[:,2]
    #get the index of the closest value of A
    index = findmin(abs.(A_data .- A))[2]
    #get the C value by interpolate between the two closest values of A
    Ai = A_data[index]
    Ai_1 = A_data[index+1] #the next value of A
    C = C_data[index] + (C_data[index+1] - C_data[index])/(Ai_1 - Ai)*(A - Ai)
    return C
end

"""
function get_Icrack input depth C output Icrack from a csv file. 
csv file dataformat is C,Icrack
"""
function get_Icrack(C::Float64)
    #read the csv file
    filename = "CtoIcrack.csv" #Make sure you got csv in the name
    data = CSV.read(filename, DataFrame)
    #get the column of C
    C_data = data[:,1]
    #get the column of Icrack
    Icrack_data = data[:,2]
    #get the index of the closest value of C
    index = findmin(abs.(C_data .- C))[2]
    #get the Icrack value
    Icrack = Icrack_data[index]
    return Icrack
end