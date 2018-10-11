cd "/Users/aaronmarkiewitz/Box/research/sra/data"
clear all
use ../data/merged_data, clear
set scheme s1color
***
* Compute regression variables
***
local beta = .99
ren fyfsgda surplus
reg surplus l.surplus, nocons
local rho = _b[l.surplus]


gen pvd_surplus = 1/(1-`beta'*`rho') * surplus

gen lol = gfdgdpa - pvd_surplus


reg surplus l.lol





