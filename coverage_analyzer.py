# Using Fountain_anallyzer.py

from Analysis.Fountain_analyzer import FT_Analyzer
from Model.Model import DNA_Channel_Model
from Model.config import DEFAULT_PASSER
import matplotlib.pyplot as plt

# # Parameters
file_name = "lena.jpg"        # your input file (put it in /files folder)
alpha = 0.3                   # coverage = 1 + alpha
rs_length = 4                 # Reed-Solomon symbols
chunk_size = 20               # data chunk size

# # Create analyzer
Model = DNA_Channel_Model(None,DEFAULT_PASSER)
analyzer = FT_Analyzer(file_name, Model, alpha=0.3, rs_length=rs_length, chunk_size=chunk_size)

# # Run encode + simulate + decode
# analyzer.run()

# # Plot decoding failure probability over different coverage levels
# analyzer.alpha_scan()
# plt.title("Decoding Failure vs. Coverage")
# plt.grid(True)
# plt.show()


# worked
for _ in range(10):
    analyzer.run()
fail_prob = analyzer.fail_prob(alpha=0.3, plot=True)
import matplotlib.pyplot as plt
plt.title("Decoding Failure at 1.3x Coverage")
plt.show()
# -