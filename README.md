****
# Project1 of FE5222
****
## About this project


- Team **HO Ngok Chao, GAO Jisan, CHENG Toyuan**

- **Summary**: This notebook in **Python 3** is to serve as part of the final report together with source code for your project, stating the logic of the derivation of pricing method (if not discussed in class).

- **Reference**: [Longstaff-Schwartz (2001): "Valuing American Options by Simulation: A Simple Least-Squares Approach." Review of Financial Studies, Vol. 14, 113-147](https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0CCEQFjAAahUKEwiXtNSZm4rHAhXHOhQKHTjBD3k&url=https%3A%2F%2Fpeople.math.ethz.ch%2F~hjfurrer%2Fteaching%2FLongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf&ei=7PO9VZeOBcf1ULiCv8gH&usg=AFQjCNFQr1r_Cf_pxylg_amU3TFOZVDc8w&sig2=ixZnX_wWQ48G66BMuQTPZA&bvm=bv.99261572,d.d24),
[Least-Square Monte Carlo for American Options in a Python Class](https://github.com/jpcolino/IPython_notebooks/blob/master/Least%20Square%20Monte%20Carlo%20Implementation%20in%20a%20Python%20Class.ipynb)

Project One
October 2, 2020
1 Requirements
1. Each group is required to submit a nal report together with source code for your project, stating in
details the derivation of pricing method (if not discussed in class), the choice of numerical algorithm,
test results and analysis of results.
2. Deadline: November 8, 2020
2 Pricing American Option
In this project, you will
1. implement either Least Square Monte Carlo (LSMC) method or Barone-Adesi and Whaley (BAW)
approximation to price American put options,
2. implement either BBS or BBSR as a benchmark model,
3. compare the results from two pricing methods you have implemented for different choices of time to
maturities, strikes and spots levels etc,
4. discuss how the choices of numerical parameters such as number of Monte Carlo sample paths in
LSMC affect the results of LSMC or BAW methods.
3 References
A useful reference for BAW approximation is the book The Complete Guide to Option Pricing Formulas
by Espen Gaarder Haug.
For LSMC method, the original paper from Longsta and Schwartz (given in the lecture notes) can
also be helpful.
