all: R/gen_data.out R/sim.Rout R/data.Rout R/examples.Rout manuscript.pdf

R/gen_data.out: R/gen_data.sh R/sim_data.R
	source gen_data.sh

R/sim.Rout: R/sim.R R/gen_data.sh R/sim_data.R
	Rscript R/sim.R 4 50 100

R/data.Rout: R/data.R data/eeesr.csv
	Rscript R/data.R 4 50 100

R/examples.Rout: R/examples.R
	Rscript R/examples.R 200

manuscript.pdf: manuscript.md
	pandoc $< -o $@ --bibliography=ref.bib -V geometry:margin=1.25in
