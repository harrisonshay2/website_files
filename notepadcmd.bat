@echo off
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\daily_update1.R"
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\daily_update1.R" 
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\daily_update1.R"
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\daily_update1.R"
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\complete_cluster_setup.R"
"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" "C:\Users\Harri\OneDrive\Documents\centauri\pure_cluster_setup.R"

"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/Country_breakdown/by_country.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\Country_breakdown\by_country.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\Country_breakdown\index.html" "C:\Users\Harri\OneDrive\Documents\Country_breakdown\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/Country_breakdown/docs"
git add .
git commit -m"by country daily"
git push

"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/Cumulative_conjunction_complete/complete_cumnum.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\Cumulative_conjunction_complete\complete_cumnum.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\Cumulative_conjunction_complete\index.html" "C:\Users\Harri\OneDrive\Documents\Cumulative_conjunction_complete\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/Cumulative_conjunction_complete/docs"
git add .
git commit -m"cumulative conjunction complete cluster daily"
git push

"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/Risk_algorithms/top50_list.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\Risk_algorithms\top50_list.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\Risk_algorithms\index.html" "C:\Users\Harri\OneDrive\Documents\Risk_algorithms\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/Risk_algorithms/docs"
git add .
git commit -m"Risk algorithms daily"
git push


"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/Conjunction_pure_cluster/pure_cumnum.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\Conjunction_pure_cluster\pure_cumnum.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\Conjunction_pure_cluster\index.html" "C:\Users\Harri\OneDrive\Documents\Conjunction_pure_cluster\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/Conjunction_pure_cluster/docs"
git add .
git commit -m"Cumulative Conjunction Pure Cluster daily"
git push

"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/pure_by_country/pure_by_country.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\pure_by_country\pure_by_country.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\pure_by_country\index.html" "C:\Users\Harri\OneDrive\Documents\pure_by_country\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/pure_by_country/docs"
git add .
git commit -m"By Country Pure Cluster daily"
git push


"C:\Program Files\R\R-4.0.0\bin\Rscript.exe" -e "library('rmarkdown');render('C:/Users/Harri/OneDrive/Documents/Risk_algorithm_method2/risk+algorithm_method2.Rmd')"
ren "C:\Users\Harri\OneDrive\Documents\Risk_algorithm_method2\risk+algorithm_method2.html" "index.html"
move "C:\Users\Harri\OneDrive\Documents\Risk_algorithm_method2\index.html" "C:\Users\Harri\OneDrive\Documents\Risk_algorithm_method2\docs\index.html"
cd "C:/Users/Harri/OneDrive/Documents/Risk_algorithm_method2/docs"
git add .
git commit -m"Statistically Most Concerning Pure Cluster daily"
git push


