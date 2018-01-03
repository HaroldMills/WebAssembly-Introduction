function timeFfts() {

    const freq = 1;
    const minExp = 4;
    const maxExp = 16;
    const maxFftCost = _fftCost(maxExp);
    const numTrials = 10;
    const minNumFftsPerTrial = 1000;

    for (let e = minExp; e <= maxExp; ++e) {

        const fftSize = 1 << e;
        const input = _getFftTestInput(fftSize, freq);
        const output = new Float64Array(fftSize);

        // Get trial size.
        const fftCostRatio = maxFftCost / _fftCost(e);
        const numFftsPerTrial = Math.round(
          minNumFftsPerTrial * fftCostRatio);

        let minUsPerFft = 1e20;

        for (let i = 0; i != numTrials; ++i) {

            const begin = window.performance.now();

            for (let j = 0; j != numFftsPerTrial; ++j)
                realFft(input, output);

            const end = window.performance.now();
            const elapsed = end - begin;
            const usPerFft = 1e3 * elapsed / numFftsPerTrial;

            if (usPerFft < minUsPerFft)
                minUsPerFft = usPerFft;

            // console.log(`    ${i} ${numFftsPerTrial} ${begin} ${end} ${elapsed} ${usPerFft}`);

        }

        // console.log(`${fftSize} ${numFftsPerTrial} ${minUsPerFft}`);
        console.log(`${minUsPerFft}`);

    }

};


function _fftCost(n) {
    return n * (1 << n);
};


function _getFftTestInput(n, f) {
    const x = new Float64Array(n);
    const factor = 2 * Math.PI * f / n;
    for (let i = 0; i != n; ++i)
        x[i] = Math.cos(factor * i);
    return x;
};
