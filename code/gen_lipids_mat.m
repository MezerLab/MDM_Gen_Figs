function lipids=gen_lipids_mat()

lipids.TotalPhos=[34.7,31.0,22.3,33.3,17.8,20.8,20.5];
lipids.TotalPhoserr=[4.3,2.6,1.9,4.7,1.2,1.8,1.9];

lipids.chol=[20.5,21.3,15.4,16.6,9.1,9.3,8.3];
lipids.cholerr=[4.7,3.1,1.1,1.4,0.9,1.2,1];

lipids.PhosTOpro=[0.315,0.31,0.223,0.333,0.228,0.214,0.225];
lipids.PhosTOproerr=[0.028,0.034,0.019,0.027,0.019,0.023,0.028];

lipids.ROI={'WM-Frontal','Pons','Cerebellum','Medulla','CTX-Frontal','Caudate Nucleus','Hippocampus'};
lipids.FS={[13,32],[18],[7,26],[19],[9,28],[2,21],[5,24]};
lipids.all={'PE [p.u.]','PS [p.u.]','PtdCho [p.u.]','PI [p.u.]','Spg [p.u.]','Phospholipids/Cholesterol [fraction]','Phospholipids/Proteins [fraction]'};

lipids.data=[45.02923977	46.92982456	44.15204678	42.69005848	34.64912281	33.33333333	35.0877193	6.578947368	7.894736842	7.456140351	5.701754386	6.140350877	8.479532164	7.894736842	26.4619883	27.63157895	31.28654971	29.97076023	44.73684211	41.37426901	37.42690058	7.01754386	6.725146199	7.01754386	8.479532164	6.432748538	7.894736842	8.771929825	15.78947368	11.40350877	10.67251462	14.03508772	8.771929825	9.795321637	11.69590643];
lipids.data=reshape(lipids.data,[7,5]);
lipids.err=[48.1865284974	50.3454231434	47.5820379965	48.3592400691	38.5146804836	36.2694300518	38.0829015544	8.0310880829	9.4127806563	8.5492227979	7.2538860104	7.6856649396	9.4127806563	9.3264248705	28.6701208981	31.2607944732	34.8013816926	32.8151986183	49.3955094991	45.5094991364	40.8462867012	8.4628670121	8.0310880829	7.9447322971	9.7582037997	7.5993091537	8.9810017271	10.3626943005	18.8255613126	12.8670120898	12.9533678756	15.8894645941	10.103626943	11.3126079447	13.0397236615];
lipids.err=reshape(lipids.err,[7,5]);
lipids.err=lipids.err-lipids.data;

lipids.percent=lipids.data;
lipids.percent(:,6)=lipids.TotalPhos./lipids.chol;
lipids.percent(:,7)=lipids.PhosTOpro;
se_x=lipids.TotalPhoserr;
se_y=lipids.cholerr;
x=lipids.TotalPhos;
y=lipids.chol;
rationerr=((y.^-2).*se_x.^2+(x./y.^2).^2.*se_y.^2);
lipids.err(:,6)=rationerr;
lipids.err(:,7)=lipids.PhosTOproerr;
end