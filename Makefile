test:
	julia --project=@. -e 'import Pkg; Pkg.test()'

coverage:
	julia --project -e 'using LocalCoverage; report_coverage_and_exit(target_coverage=90)'

clean:
	julia --project=@. -e 'using LocalCoverage; clean_coverage()'

build:
	julia --project=@. -e 'import Pkg; Pkg.build()'

.PHONY: coverage clean build test
