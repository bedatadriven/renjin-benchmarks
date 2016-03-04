$R -q -e 'install.packages("devtools")'
$R -q -e 'devtools::install_github("rstudio/packrat")'
$R -q -e 'source("http://bioconductor.org/biocLite.R")'
# loop through directories
for f in $( find -path './[^.]*' -prune -type d ); do
  echo "-----------------------------------------------"
  # print directory name
  echo $f
  #set counter to 0
  COUNT=0
    while [ $COUNT -le $1 ];
      do
          count_info=`printf "round %s out of %s | warm runs: %s " "$COUNT" "$1" "$2"`
          echo $count_info
          cd $f
          $R -q -e "print(getwd());"
          $R -q -e 'library(packrat); packrat::restore();'
          $R -q -e "library(benchmarkR); runBenchmark($2); db <- Sys.getenv(c('BENCH_DB', 'BENCH_TYPE')); benchDBReport(db_name = db[[1]], con_type = db[[2]] )"
          cd ..
          COUNT=$(( $COUNT + 1 ))
      done
  echo "-----------------------------------------------"
done