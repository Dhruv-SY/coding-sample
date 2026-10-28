****************************************************
* 17_graphs.do
* RD graphs in the style of Anderson (2018), Figures 2-4:
* outcome means in 5km distance bins over +/-100km, a
* linear fit on each side of the border, a red vertical
* line at the border, and side labels inside the plot.
****************************************************

clear all
set more off

*---------------------------------------------------*
* 0) CALL SHARED SETUP
*---------------------------------------------------*
do "code/00_setup.do"

*---------------------------------------------------*
* 1) GRAPH SETUP
*    Fixed-width 5km bins over +/-100km (Anderson 2018 uses
*    5km bins over the same range); bin midpoints for plotting
*---------------------------------------------------*

capture drop sample100 bin_mid
gen sample100 = rd > -100 & rd < 100
gen bin_mid = floor(rd/5)*5 + 2.5 if sample100

*---------------------------------------------------*
* 2) REUSABLE GRAPH PROGRAM
*    NOTE: these pool all 12 border pairs without pair fixed
*    effects, so bin means partly reflect PAIR COMPOSITION,
*    not just the discontinuity. The regressions (which include
*    border-pair FE and pair-specific slopes) are immune to
*    this; treat these figures as descriptive, not as visual
*    estimates of the RD coefficient.
*---------------------------------------------------*

capture program drop rd_graph
program define rd_graph
    syntax varname(numeric), YTITLE(string) FILE(string)

    capture drop ymean_graph tag_graph

    egen ymean_graph = mean(`varlist') if sample100 & !missing(`varlist'), by(bin_mid)
    egen tag_graph   = tag(bin_mid)    if sample100 & !missing(`varlist')

    * y position for the side labels: just above the highest bin mean,
    * with the y-scale extended so the labels stay inside the plot
    quietly summarize ymean_graph if tag_graph == 1
    local yrange = r(max) - r(min)
    local ytop   = r(max) + 0.08*`yrange'
    local yhigh  = r(max) + 0.14*`yrange'
    local ylow   = r(min)

    twoway ///
        (scatter ymean_graph bin_mid if tag_graph == 1, ///
            mcolor(navy) msymbol(diamond) msize(small)) ///
        (lfit `varlist' rd if rd < 0 & rd > -100 & !missing(`varlist'), ///
            lcolor(navy) lwidth(medium)) ///
        (lfit `varlist' rd if rd >= 0 & rd < 100 & !missing(`varlist'), ///
            lcolor(navy) lwidth(medium)) ///
        , ///
        xline(0, lcolor(red) lwidth(medthin)) ///
        text(`ytop' -50 "Laxer law", size(small) color(gs5)) ///
        text(`ytop'  50 "Stricter law", size(small) color(gs5)) ///
        xtitle("Distance to border") ///
        ytitle("`ytitle'") ///
        yscale(range(`ylow' `yhigh')) ///
        xlabel(-100(50)100, nogrid) ///
        ylabel(#5, angle(horizontal) glcolor(gs14) glwidth(thin)) ///
        legend(off) ///
        graphregion(color(white)) ///
        plotregion(color(white))

    graph export "results/`file'", replace width(2000)

    capture drop ymean_graph tag_graph
end

****************************************************
* 3) RD GRAPHS
****************************************************

rd_graph child_marriage, ///
    ytitle("Child marriage") ///
    file("child_marriage_rd_plot.png")

rd_graph multiple_unions, ///
    ytitle("Multiple unions") ///
    file("multiple_unions_rd_plot.png")

rd_graph v133, ///
    ytitle("Education (years)") ///
    file("education_years_rd_plot.png")

rd_graph educ_any, ///
    ytitle("Any education") ///
    file("any_education_rd_plot.png")

rd_graph educ_years_pos, ///
    ytitle("Education years among schooled women") ///
    file("education_years_positive_rd_plot.png")

rd_graph v106, ///
    ytitle("Education level") ///
    file("education_level_rd_plot.png")

rd_graph v212, ///
    ytitle("Age at first birth") ///
    file("age_first_birth_rd_plot.png")

rd_graph v201, ///
    ytitle("Children ever born") ///
    file("children_ever_born_rd_plot.png")

rd_graph v525, ///
    ytitle("Age at first sex") ///
    file("age_first_sex_rd_plot.png")

rd_graph not_working, ///
    ytitle("Not working") ///
    file("not_working_rd_plot.png")

rd_graph domestic_work, ///
    ytitle("Domestic work among workers") ///
    file("domestic_work_rd_plot.png")

rd_graph paid_cash, ///
    ytitle("Paid in cash among workers") ///
    file("paid_cash_rd_plot.png")

rd_graph white_collar, ///
    ytitle("White-collar employment") ///
    file("white_collar_rd_plot.png")

rd_graph white_collar_cond, ///
    ytitle("White-collar employment among workers") ///
    file("white_collar_conditional_rd_plot.png")
