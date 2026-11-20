import Link from 'next/link';
import earlier from '@/public/data/omnimodel-reports/k/2796-kaprox-mod.json';
import later from '@/public/data/omnimodel-reports/k/101629-kaprox-mod.json';
import { comparisonPath } from '@/lib/comparison';
import styles from './home.module.css';

const example = comparisonPath('2796_kaprox-ID-1', '101629_kaprox-ID-1');
const models = [
  { id: '2796', year: 1999, report: earlier, color: '#9ce6bb' },
  { id: '101629', year: 2008, report: later, color: '#ffb781' },
];

function activationPoints({ id, report }: (typeof models)[number]) {
  if (report.channel !== 'K' || report.originalSlug !== id + '_kaprox.mod') {
    throw new Error('Unexpected model in kaprox comparison');
  }
  function series(key: string): number[] {
    const value = report.sections.find(section => section.key === key)?.value;
    if (value && typeof value === 'object' && 'n' in value && Array.isArray(value.n)
      && value.n.length > 0 && value.n.every(Number.isFinite)) return value.n;
    throw new Error('Missing kaprox activation series: ' + key);
  }
  const voltage = series('RATE_VALS_V');
  const activation = series('RATE_VALS_SS');
  if (voltage.length !== activation.length) throw new Error('Mismatched kaprox activation samples');
  return voltage.map((v, i) => [x(v), y(activation[i])].join(',')).join(' ');
}

const x = (mv: number) => 44 + (mv + 100) * 2.08;
const y = (fraction: number) => 187 - fraction * 166;

function ActivationPlot() {
  return <figure className={styles.figure}>
    <figcaption><span>K<sup>+</sup> &nbsp;/&nbsp; kaprox.mod</span><strong>Steady-state activation</strong></figcaption>
    <svg viewBox="0 0 485 244" role="img" aria-label="Steady-state activation of gate n against membrane voltage. The activation curve for ModelDB 101629 lies to the right of the curve for 2796.">
      {[0, 0.25, 0.5, 0.75, 1].map(v => <g key={v}>
        <line x1="44" x2="460" y1={y(v)} y2={y(v)} stroke="#426258" strokeDasharray="2 5" />
        <text x="33" y={y(v) + 4} textAnchor="end">{v.toFixed(2)}</text>
      </g>)}
      {[-100, -50, 0, 50, 100].map(v => <g key={v}>
        <line x1={x(v)} x2={x(v)} y1="21" y2="187" stroke="#38574e" strokeDasharray="2 5"/>
        <text x={x(v)} y="209" textAnchor="middle">{v}</text>
      </g>)}
      {models.map(m => <polyline key={m.id} points={activationPoints(m)} fill="none" stroke={m.color}
        strokeWidth="3" strokeDasharray={m.id === '101629' ? '7 4' : undefined} strokeLinejoin="round"/>)}
      <text x="252" y="237" textAnchor="middle">Membrane voltage (mV)</text>
    </svg>
    <div className={styles.legend}>{models.map(m =>
      <a key={m.id} href={'/omnimodels/K/' + m.report.originalSlug} style={{color:m.color}}>
        <i style={{borderTopStyle:m.id === '101629' ? 'dashed' : 'solid'}} />ModelDB {m.id} · {m.year}
      </a>)}</div>
  </figure>;
}

export default function Home() {
  return <div className={styles.page}><div className={styles.wrap}>
    <header className={styles.mast}><Link href="/" className={styles.brand}>ICG <em>Explorer</em></Link>
      <nav aria-label="Main navigation"><Link href="/visualizer">Graph</Link><Link href="/omnimodels">Omnimodels</Link></nav>
    </header>
    <main>
      <section className={styles.hero}>
        <div>
          <h1>Ion Channel<br/><em>Genealogy</em></h1>
          <p className={styles.intro}>Trace how researchers reuse and modify ion-channel models across ModelDB.</p>
          <p className={styles.case}>Two potassium-channel files share <strong>98.75% code similarity</strong>. The later version adds a voltage shift to the gating equations. See where the source changes and how the activation curves differ.</p>
          <Link className={styles.cta} href={example}>Compare kaprox.mod <span aria-hidden="true">→</span></Link>
        </div>
        <ActivationPlot />
      </section>
      <section className={styles.paths} aria-label="Explore the collection">
        <Link href="/visualizer"><span className={styles.count}>4,820 models</span><h2>Explore the graph <span aria-hidden="true">↗</span></h2><p>Follow 38,100 source-similarity links across the collection.</p></Link>
        <Link href="/omnimodels"><span className={styles.count}>3,524 reports</span><h2>Browse Omnimodels <span aria-hidden="true">↗</span></h2><p>Inspect fitted parameters, steady states and time constants.</p></Link>
      </section>
    </main>
    <footer className={styles.footer}><span>Scientific sources</span><a href="https://modeldb.science">ModelDB ↗</a><a href="https://icg.neurotheory.ox.ac.uk/">Ion Channel Genealogy ↗</a></footer>
  </div></div>;
}
