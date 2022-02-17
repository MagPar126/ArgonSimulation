from MolDyn import MolecularDynamics

model = MolecularDynamics(1, 2, 0.01, dimension=3)

print(model.minimum_image_convention())


