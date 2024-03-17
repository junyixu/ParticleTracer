# y = [[GenRandNum_Para(1.0), 2π*rand(), 2π*rand()] for i in 1:1000]
# yy = y .|> x->TORD2CART(x, UserInputs.R0)
# new=stack(yy)
# x = new[1, :]
# y = new[2, :]
# z = new[3, :]
# R = @. sqrt(x^2 + y^2)
# plt.scatter(R, z)

R0=UserInputs.R0

# %%
arr = [SetParticlePosition_ParabolicTorus(0.4/u.x) for i = 1:10000]
R_(arr) = sqrt(arr[1]^2 + arr[2]^2)
Z_(arr) = arr[3]
R = R_.(arr)
Z = Z_.(arr)
plt.scatter(R .- 1.7, Z)

XX = stack(arr)

ax=plt.axes(projection="3d")
x = XX[1, :]
y = XX[2, :]
z = XX[3, :]
ax.scatter3D(x, y, z)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

# %%
XX = XX / u.x
PP = [SetParticleMomentum_Gyrocenter(XX[:, i]..., Fields.tokamak(XX[:, i]..., 1.0)) for i = 1:length(XX[1, :])]
E_.(PP)
mean(norm.(PP))
