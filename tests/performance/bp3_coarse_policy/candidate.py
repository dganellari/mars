"""Prepare the isolated larger-root experiment from the current production sources."""
from pathlib import Path


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Source changed at: ' + old)
    return text.replace(old, new)


def generate(source, out):
    out.mkdir(parents=True, exist_ok=False)
    adapter = (source / 'backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu').read_text()
    # Keep adaptive -1 identical. Only explicit fixed limits can reach the larger level.
    for old, new in (
        ('requested<=50000', 'requested<=500000'),
        ('result<=50000', 'result<=500000'),
        ('limit<=50000', 'limit<=500000')):
        count = 2 if old == 'limit<=50000' else 1
        if adapter.count(old) != count:
            raise ValueError('Unexpected cutoff guards')
        adapter = adapter.replace(old, new)
    adapter = adapter.replace('between 10 and 50000 rows', 'between 10 and 500000 rows')
    needle = '        ++audit.setups; audit.empty_setups+=(local==0);'
    adapter = replace_once(adapter, needle, '''#ifdef BP3_CONSOLIDATE_CUDA
        if (rank==0) {
            size_t free_bytes=0,total_bytes=0;
            gpu(cudaMemGetInfo(&free_bytes,&total_bytes));
            std::printf("BP3 root memory: %zu used of %zu bytes after setup\\n",total_bytes-free_bytes,total_bytes);
        }
#endif
''' + needle)
    folder = out / 'solvers'
    folder.mkdir()
    (folder / 'mars_cellwise_amg_consolidate.cu').write_text(adapter)
    driver = (source / 'examples/distributed/unstructured/mars_marsir_bp.cu').read_text()
    driver = replace_once(driver, 'rows < -1 || rows > 50000', 'rows < -1 || rows > 500000')
    driver = replace_once(driver, 'o.consolidate_rows > 50000', 'o.consolidate_rows > 500000')
    driver = driver.replace('10..50000', '10..500000')
    folder = out / 'driver'
    folder.mkdir()
    (folder / 'mars_marsir_bp.cu').write_text(driver)
    fixture = (source / 'tests/bp3_coarse/fixture.cu').read_text()
    fixture = replace_once(fixture, 'backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu',
                           'solvers/mars_cellwise_amg_consolidate.cu')
    fixture = replace_once(fixture, 'n<=32767', 'n<=262143')
    (out / 'fixture.cu').write_text(fixture)
