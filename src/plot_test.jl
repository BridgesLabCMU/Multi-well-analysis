using Plots
using LaTeXStrings

pgfplotsx()

push!(PGFPlotsX.CUSTOM_PREAMBLE,
	  """
\\usepackage[scaled]{helvet}
\\renewcommand\\familydefault{\\sfdefault}
\\usepackage[T1]{fontenc}
\\usepackage{helvet, sansmath}
\\sansmath
""")

x = range(0,10,length=100)
y = sin.(x)
y_noisy = @. sin(x) + 0.1*randn()

plot(x, y, label=latexstring("sin(x)"), lw=2)
scatter!(x, y_noisy, label="data", mc=:red, ms=2, ma=0.5)
plot!(legend=:bottomleft)
title!("Sine with noise, plotted with GR")
xlabel!("x")
ylabel!("y")
savefig("test.pdf")
