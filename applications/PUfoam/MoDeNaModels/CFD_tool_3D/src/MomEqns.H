if
(
    PBESwitch
 && PBEMethod == "QMOM"
 && nNodes == 2
)
{
    fvScalarMatrix mZeroEqn
    (
        fvm::ddt(mZero)
      + fvm::div(phi,mZero)
      + fvm::SuSp(expImpSource,mZero)
      - fvm::div(alpha1phir,mZero)
     == mZeroSource
    );
    mZeroEqn.relax();

    fvScalarMatrix mOneEqn
    (
        fvm::ddt(mOne)
      + fvm::div(phi,mOne)
      + fvm::SuSp(expImpSource,mOne)
      - fvm::div(alpha1phir,mOne)
     == mOneSource
    );
    mOneEqn.relax();

    fvScalarMatrix mTwoEqn
    (
        fvm::ddt(mTwo)
      + fvm::div(phi,mTwo)
      + fvm::SuSp(expImpSource,mTwo)
      - fvm::div(alpha1phir,mTwo)
     == mTwoSource
    );
    mTwoEqn.relax();

    fvScalarMatrix mThreeEqn
    (
        fvm::ddt(mThree)
      + fvm::div(phi,mThree)
      + fvm::SuSp(expImpSource,mThree)
      - fvm::div(alpha1phir,mThree)
     == mThreeSource
    );
    mThreeEqn.relax();

    Info<< "\nSolving Moments Equations. " << endl;
    mZeroEqn.solve();
    mOneEqn.solve();
    mTwoEqn.solve();
    mThreeEqn.solve();
}

if
(
    PBESwitch
 && PBEMethod == "QMOM"
 && nNodes == 3
)
{
    fvScalarMatrix mZeroEqn
    (
        fvm::ddt(mZero)
      + fvm::div(phi,mZero)
      + fvm::SuSp(expImpSource,mZero)
      - fvm::div(alpha1phir,mZero)
     == mZeroSource
    );
    mZeroEqn.relax();

    fvScalarMatrix mOneEqn
    (
        fvm::ddt(mOne)
      + fvm::div(phi,mOne)
      + fvm::SuSp(expImpSource,mOne)
      - fvm::div(alpha1phir,mOne)
     == mOneSource
    );
    mOneEqn.relax();

    fvScalarMatrix mTwoEqn
    (
        fvm::ddt(mTwo)
      + fvm::div(phi,mTwo)
      + fvm::SuSp(expImpSource,mTwo)
      - fvm::div(alpha1phir,mTwo)
     == mTwoSource
    );
    mTwoEqn.relax();

    fvScalarMatrix mThreeEqn
    (
        fvm::ddt(mThree)
      + fvm::div(phi,mThree)
      + fvm::SuSp(expImpSource,mThree)
      - fvm::div(alpha1phir,mThree)
     == mThreeSource
    );
    mThreeEqn.relax();

    fvScalarMatrix mFourEqn
    (
        fvm::ddt(mFour)
      + fvm::div(phi,mFour)
      + fvm::SuSp(expImpSource,mFour)
      - fvm::div(alpha1phir,mFour)
     == mFourSource
    );
    mFourEqn.relax();

    fvScalarMatrix mFiveEqn
    (
        fvm::ddt(mFive)
      + fvm::div(phi,mFive)
      + fvm::SuSp(expImpSource,mFive)
      - fvm::div(alpha1phir,mFive)
     == mFiveSource
    );
    mFiveEqn.relax();

    Info<< "\nSolving Moments Equations. " << endl;
    mZeroEqn.solve();
    mOneEqn.solve();
    mTwoEqn.solve();
    mThreeEqn.solve();
    mFourEqn.solve();
    mFiveEqn.solve();
}