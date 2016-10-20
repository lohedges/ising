#!/usr/bin/awk -f
# ising.awk: a simple 2D Ising model simulation in awk

# Copyright (C) 2013  Lester Hedges <lester.hedges+ising@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ising.awk supports the following command-line arguments, defaults
# are given in the BEGIN action.
#
# t=, temperature=
#           Temperature of the thermal bath.
# w=, width=
#           Width of the lattice.
# j=, J=
#           Ising coupling constant.
# h=, H=
#           External field.
# s=, sweeps=
#           Number of Monte Carlo sweeps.
# i=, initialize=
#           Initial condition, random, spin down, or spin up.
# o=, output=
#           Frequency of output in Monte Carlo sweeps.
# r=, restart=
#           Restart file (to be loaded).
# c=, config=
#           Restart file name (for output).
# l=, logfile=
#           File name for logging.
# wolff=
#           Whether to use Wolff algorithm ("on" or "off")
#           See http://en.wikipedia.org/wiki/Wolff_algorithm

# Since there are no files to process everything is done within
# the BEGIN action. The order of operation is as follows:
#   1) initialize defaults
#   2) parse any command-line arguments
#   3) execute main Ising model simulation
BEGIN {
    # set defaults
    width = 60
    sites = width*width
    temperature = 2.269
    beta = 1/temperature
    J = 1
    h = 0
    max_sweeps = 1e4
    initialize = "random"
    report_interval = 100*sites
    log_file = "report.txt"
    restart_file = "restart.txt"
    is_wolff = 0

    # seed random number generator
    srand()

    # parse arguments from command-line
    parse_cmd_line_args()

    # execute main Ising simulation
    main()
}

# Parse command-line arguments
function parse_cmd_line_args(  i, s) {
    for (i=1;i<ARGC;i++) {
        split(ARGV[i],s,"=")

        if (s[1] == "t" || s[1] == "temperature") {
            temperature = s[2]
            beta = 1/temperature
        }
        else if (s[1] == "w" || s[1] == "width") {
            width = s[2]
            sites = width*width
            sites_set = 1
        }
        else if (s[1] == "j" || s[1] == "J") {
            J = s[2]
        }
        else if (s[1] == "h" || s[1] == "H") {
            h = s[2]
        }
        else if (s[1] == "s" || s[1] == "sweeps") {
            max_sweeps = s[2]
        }
        else if (s[1] == "i" || s[1] == "initialize") {
            initialize = s[2]
        }
        else if (s[1] == "o" || s[1] == "output") {
            report_interval = s[2]
            report_interval_set = 1
        }
        else if (s[1] == "r" || s[1] == "restart") {
            load_lattice(s[2])
            restart_set = 1
        }
        else if (s[1] == "c" || s[1] == "config") {
            restart_file = s[2]
        }
        else if (s[1] == "l" || s[1] == "logfile") {
            log_file = s[2]
        }
        else if ("wolff") {
            if (s[2] == "on") is_wolff=1
        }
        else {
            print "Error: unknown command line argument \""s[1]"\"."
            exit 1
        }
    }

    if (restart_set != 1) {
        if (sites_set = 1) {
            if (report_interval_set == 1) {
                report_interval *= sites
            }
            else report_interval = 100*sites
        }
        else report_interval *= sites
    }
    else {
        if (sites_set == 1) {
            if (restart_width != width) {
                printf("Warning: lattice width mismatch, ")
                printf("using restart configuration size.\n")
                width = restart_width
                sites = width*width
            }
            if (report_interval_set) report_interval *= sites
            else report_interval = 100*sites
        }
        else {
            width = restart_width
            sites = width*width

            if (report_interval_set) report_interval *= sites
            else report_interval = 100*sites
        }
    }

    # set link probability for is_wolff algorithm
    if (is_wolff == 1) link_prob = 1.0 - exp(-2*beta)
}

# Main function: execute Ising simulation
function main() {
    # generate random starting configuration
    if (restart_set == 0) {
        initialize_lattice()
    }

    # compute initial energy
    total_energy = get_total_energy()

    # incremental time step (in MC sweeps)
    t_step = 1/sites

    # rescale time for cluster moves
    if (is_wolff) {
        report_interval /= sites
        t_step = 1
    }

    # print simulation info (parameter values, etc.)
    printf("Starting Ising simulation...\n")
    printf("width = %d\n", width)
    printf("temperature = %5.4f\n", temperature)
    printf("J = %5.4f\n", J)
    printf("h = %5.4f\n", h)
    printf("sweeps = %5.4e\n", max_sweeps)
    if (restart_set != 1) {
        printf("initialize = %s\n", initialize)
    }
    printf("logfile = %s\n", log_file)
    printf("images = ./Config_*.ps\n")

    printf("\nReporting every %d sweeps.\n", report_interval/sites)
    if (is_wolff) print "Using Wolff algorithm."

    # initialize output files and formatting
    initialize_output()

    # MAIN LOOP
    while (sweeps < max_sweeps) {
        sweeps += t_step
        report_flag++

        if (is_wolff == 0) metropolis_step()
        else wolff()

        # write report log to file and stdout then generate
        # a postscript of the current lattice configuration
        if (report_flag == report_interval) {
            # recompute energy and magnetization if using Wolff algorithm
            if (is_wolff) total_energy = get_total_energy()
            frame++
            report(frame)
            report_flag = 0
        }
    }
    # END MAIN LOOP

    exit 0
}

# Initialize spin configuration
function initialize_lattice(  i) {
    magnetization = 0
    for (i=0;i<sites;i++) {
        if (initialize == "random") {
            rand() < 0.5 ? lattice[i] = -1 : lattice[i] = 1
        }
        else if (initialize == "up" || initialize == "+") {
            lattice[i] = 1
        }
        else if (initialize == "down" || initialize == "-") {
            lattice[i] = -1
        }
        else {
            printf("Error: unknown initialization option, allowed values ")
            printf("are \"random\", \"up\", and \"down\"\n")
            exit 1
        }
        magnetization += lattice[i]
    }
}

# Load spin configuration from file
function load_lattice(f,  x, s, i) {
    tmp = -1

    while (getline x < f > 0) {

        # compute width of line
        restart_width = length(x)

        if (tmp > 0) {
            if (restart_width != tmp) {
                print "Error: line length mismatch in configuration file!"
                exit 1
            }
        }

        split(x, s, "")

        for (i=1;i<=restart_width;i++) {
            s[i] == "+" ? lattice[i - 1 + row*restart_width] = 1 :
                lattice[i -1 + row*restart_width] = -1
        }

        row++
    }

    # check that lattice is square
    if (row != restart_width) {
        print "Error: lattice isn't square. Aborting!"
        exit 1
    }
    else if (row == 0) {
        print "Error: configuration file is empty. Aborting!"
        exit 1
    }
}

# Print current spin configuration to file
function save_lattice(  i) {
    for (i=0;i<sites;i++) {
        if (i == 0) {
            printf("%s%s", lattice[i] == 1 ? "+" : "-",
                   (i+1)%width == 0 ? "\n" : "") > restart_file
        }
        else printf("%s%s", lattice[i] == 1 ? "+" : "-",
             (i+1)%width == 0 ? "\n" : "") >> restart_file
    }

    close(restart_file)
}

# Return one of the four neighbors of lattice site "s"
function get_neighbor(s, n,  x, y) {
    y = int(s/width)
    x = s - width*y

    if (n == 0) return (x - 1 + width)%width + y*width
    if (n == 1) return (x + 1)%width + y*width
    if (n == 2) return x + width*((y - 1 + width)%width)
    if (n == 3) return x + width*((y + 1)%width)
}

# Compute and return the energy of lattice site "s"
function get_energy(s,  i, e) {
    e -= h*lattice[s]
    for (i=0;i<4;i++) {
        e -= J*lattice[s]*lattice[get_neighbor(s,i)]
    }
    return e
}

# Compute and return the total lattice energy
function get_total_energy(  i, j, e) {
    # reset magnetization
    magnetization = 0

    for (i=0;i<sites;i++) {
        magnetization += lattice[i]
        e -= h*lattice[i]
        for (j=0;j<4;j++) {
            e -= J*lattice[i]*lattice[get_neighbor(i,j)]
        }
    }
    return 0.5*e
}

# Perform a Metropolis Monte Carlo trial spin flip
function metropolis_step() {
    # choose random lattice site
    site = int(sites*rand())

    # compute initial site energy
    initial_energy = get_energy(site)

    # flip site
    lattice[site] *= -1

    # compute final energy
    final_energy = get_energy(site)

    # compute energy change
    energy_change = final_energy - initial_energy

    # reset acceptance flag
    is_accepted = 0

    # check if move is accepted
    if (energy_change > 0) {
        if (rand() < exp(-beta*energy_change))
            is_accepted = 1
    }
    else is_accepted = 1

    # update energy and magnetization,
    # reset spin if trial move failed
    if (is_accepted == 1) {
        total_energy += energy_change
        magnetization += 2*lattice[site]
        accepts++
    }
    else lattice[site] *= -1
}

# Perform iteration of Wolff algorithm
function wolff() {
    # choose random lattice site
    site = int(sites*rand())

    # store spin value
    spin = lattice[site]

    cluster_size = 0

    # recursively flip cluster of spins
    recursive_cluster_flip(site, spin)
}

# Recursively flip cluster of spins
function recursive_cluster_flip(site, spin,  i, n) {
    cluster_size++

    # flip spin
    lattice[site] *= -1

    # check all neighbors of site
    for (i=0;i<4;i++)
    {
        n = get_neighbor(site, i)

        # check spin is in the same state
        if (lattice[n] == spin) {
            # add neigbor to cluster with linking probability
            if (rand() < link_prob) {
                recursive_cluster_flip(n, spin)
            }
        }
    }
}

# Wipe existing log file and write headers to file and stdout
function initialize_output() {
    printf("#------------------------------------\n") > log_file
    printf("#%9s %9s %16s\n", "sweeps", "energy", "magnetization") >> log_file
    printf("#------------------------------------\n") >> log_file
    printf("\n#------------------------------------\n")
    printf("#%9s %9s %16s\n", "sweeps", "energy", "magnetization")
    printf("#------------------------------------\n")
}

# Write energy and magnetization to stdout/file and generate
# a postscript image of the current lattice configuration
function report(f) {
    printf("%9.4e %9.4f %16.4f\n",
        sweeps, total_energy/sites, magnetization/sites)
    printf("%9.4e %9.4f %16.4f\n",
        sweeps, total_energy/sites, magnetization/sites) >> log_file
    close(log_file)
    plot_lattice(f)
    save_lattice()
}

# Write current spin configuration as a postscript
function plot_lattice(f,  i) {
    # work out scale factor
    scale = 600/width

    # generate file name
    file_name = sprintf("Config_%04d.ps", f)

    # write postscript header to file
    print "%!PS-Adobe-2.0 EPSF-2.0" > file_name
    print "%%BoundingBox: 0 0 600 600" >> file_name

    # write postscript primitive definitions
    print "/tr{translate}def" >> file_name
    print "/s {" >> file_name
    print "newpath" >> file_name
    print "0 0 moveto" >> file_name
    print "0 1 rlineto" >> file_name
    print "1 0 rlineto" >> file_name
    print "0 -1 rlineto" >> file_name
    print "closepath" >> file_name
    print "1 0 0 setrgbcolor" >> file_name
    print "fill" >> file_name
    print "} def" >> file_name

    # print scale factor to file
    printf("%5.4f %5.4f scale\n", scale, scale) >> file_name

    # set initial coordinates
    xx = 0
    yy = 0

    for (i=0;i<sites;i++) {
        if (lattice[i] == 1) {
            y = int(i/width)
            x = i - width*y

            printf("%.2f %.2f tr s\n", x-xx, y-yy) >> file_name

            # update coordinates
            xx = x
            yy = y
        }
    }

    print "showpage" >> file_name
    close(file_name)

    # copy current config to Config.ps,
    # this allows you to watch the time evolution of the lattice
    # by leaving Config.ps open in your postscript viewer
    system_cmd = sprintf("cp %s Config.ps", file_name)
    system(system_cmd)
}
