test:
	julia --project=@. -e 'import Pkg; Pkg.test()'

coverage:
	julia --project -e 'using LocalCoverage; report_coverage_and_exit(target_coverage=90)'

xcov:
	julia --project -e 'using LocalCoverage; html_coverage(open=true,dir="coverage")'

clean:
	julia --project=@. -e 'using LocalCoverage; clean_coverage()'

build:
	julia --project=@. -e 'import Pkg; Pkg.build()'

docs:
	julia --project=@. docs/make.jl

session: build
	julia --project=@.

.PHONY: test coverage clean build session docs
