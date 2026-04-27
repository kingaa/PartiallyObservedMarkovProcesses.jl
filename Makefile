JULIA = julia -t auto --project

default: build

test: build
	$(JULIA) -e 'import Pkg; Pkg.test()'

coverage: build
	$(JULIA) -e 'using LocalCoverage; report_coverage_and_exit(target_coverage=90)'

vcov: build
	$(JULIA) -e 'using LocalCoverage; html_coverage(open=true,dir="coverage")'

clean:
	$(JULIA) -e 'using LocalCoverage; clean_coverage()'
	$(RM) profile.pb.gz

build:
	$(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.build(); Pkg.precompile()'
	cd test; $(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
	cd docs; $(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'

resolve:
	$(JULIA) -e 'import Pkg; Pkg.resolve()'
	cd test; $(JULIA) -e 'import Pkg; Pkg.resolve()'
	cd docs; $(JULIA) -e 'import Pkg; Pkg.resolve()'

update:
	$(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.update(); Pkg.gc()'
	cd test; $(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.update(); Pkg.gc()'
	cd docs; $(JULIA) -e 'import Pkg; Pkg.instantiate(); Pkg.update(); Pkg.gc()'

docs: build
	make -C docs

serve: build
	$(JULIA) -e 'import LiveServer; LiveServer.servedocs()'

session: build
	$(JULIA)

opto: build
	$(JULIA) -O3

.PHONY: test coverage clean build session docs update
