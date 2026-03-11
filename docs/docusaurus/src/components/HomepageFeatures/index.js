import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Modern C++',
    Svg: require('@site/static/img/undraw_docusaurus_mountain.svg').default,
    description: (
      <>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
        Proin ut metus sit amet arcu vehicula egestas. 
        Quisque nunc elit, sagittis sed vehicula eget, egestas non eros. 
      </>
    ),
  },
  {
    title: 'AMR',
    Svg: require('@site/static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>
        Aenean ante ligula, tincidunt sed elit eu, fermentum convallis enim. 
        Duis tincidunt dui vitae arcu ultricies, sed accumsan dolor ornare.
      </>
    ),
  },
  {
    title: 'MPI',
    Svg: require('@site/static/img/undraw_docusaurus_react.svg').default,
    description: (
      <>
        Quisque diam lectus, semper ac interdum vitae, hendrerit id metus. 
        Nam enim ligula, varius in mauris ut, eleifend commodo risus. 
        Phasellus et nibh quis risus consectetur vehicula.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} role="img" />
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
