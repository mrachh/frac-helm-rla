function [cvals] = struve101(z)

    sz = size(z);
    z = z(:);
    iloc = find(abs(z)<95);
    ifar = find(abs(z)>=95);
    cvals = zeros(size(z));
    cvals(iloc) = eval_struve_loc(z(iloc));
    cvals(ifar) = eval_struve_asymp(z(ifar));
    
    cvals = reshape(cvals,sz);

end

function [cvals] = eval_struve_loc(z)

        if (numel(z) == 0)
            cvals = zeros(size(z));
            return
        end

        roots = [
       0.23287797554979764574239513845863310D-02,
       0.77797210528565602202196751376561240D-02,
       0.16271618033986618198702811217156061D-01,
       0.27707503723264139261127948179293856D-01,
       0.41963012934216623038977890242892569D-01,
       0.58891418218743405778560394085441146D-01,
       0.78329391689422605486978192014671943D-01,
       0.10010263860402900542980504967723438D+00,
       0.12403097627111796188202638019389090D+00,
       0.14993260442201184905565140071361611D+00,
       0.17762745410598842069997492381278595D+00,
       0.20693961217903811539183343424801665D+00,
       0.23769889224165259793060177883684293D+00,
       0.26974166243662189627173897294466337D+00,
       0.30291105282421662136522246065862774D+00,
       0.33705665857467929257135214342578155D+00,
       0.37203383756330568046403009390787936D+00,
       0.40770267783662548059749513310209640D+00,
       0.44392668553654047387383102654870569D+00,
       0.48057121923925020977835009680085893D+00,
       0.51750167311684907657185822394561155D+00,
       0.55458138900459695607645970207925428D+00,
       0.59166925629421969496594497355751397D+00,
       0.62861693883643152084181844298409790D+00,
       0.66526565104397613812338690028125968D+00,
       0.70144239471352053795506975786674479D+00,
       0.73695557165015248739431902109333251D+00,
       0.77158992088379493410487935842608357D+00,
       0.80510082300500561137518182658171708D+00,
       0.83720822069456348292247823421660014D+00,
       0.86759080958908695616059731434000299D+00,
       0.89588187731948532410575618640971064D+00,
       0.92166932590546832099280447416876834D+00,
       0.94450394921279951711390249056169951D+00,
       0.96392127329142721670042025053186995D+00,
       0.97948100688506504715664233260675227D+00,
       0.99082007648743246982849951429630583D+00,
       0.99769783994110755421389441552630474D+00];

        weights = [
       0.39054246308472108537572143018542489D-02,
       0.69851156794182489718899661398675367D-02,
       0.99824773011900623033011655549735733D-02,
       0.12868477385036746802691008127225257D-01,
       0.15617988869632857012575714934677317D-01,
       0.18211505504860765237462396137799908D-01,
       0.20635307068313352195120159032980171D-01,
       0.22881084883173884340069747783455705D-01,
       0.24945234599334783449814705461144249D-01,
       0.26827969985166285752458832247831409D-01,
       0.28532384616938756151580122811620477D-01,
       0.30063553632169448782541566530842284D-01,
       0.31427731341652938389422977827806917D-01,
       0.32631669241865390472468304928884732D-01,
       0.33682055948710753297697772662680929D-01,
       0.34585065860587188374649279101325372D-01,
       0.35345995387008851284279099044910354D-01,
       0.35968962282864375803451612771640899D-01,
       0.36456643112443903367680001996253688D-01,
       0.36810024688899325062925443352680342D-01,
       0.37028146561546252296714166666237650D-01,
       0.37107812810054908793043848494316541D-01,
       0.37043252593480939850414399158839917D-01,
       0.36825710732858259789984880122231590D-01,
       0.36442953627131857006501666929733034D-01,
       0.35878685118130530049529049930153823D-01,
       0.35111887371695664011105788319219062D-01,
       0.34116143950671667563908738183976545D-01,
       0.32859083990848582685090102345858933D-01,
       0.31302235858241964317952738389628315D-01,
       0.29401832538196193841361929053862802D-01,
       0.27111495402561458561401180300153319D-01,
       0.24388182350537953056835591739052492D-01,
       0.21202989623306897961392790227859606D-01,
       0.17557342613948635157400627184362894D-01,
       0.13500760099128853258969852643313634D-01,
       0.91371732847652890752352921528928191D-02,
       0.46001586784217374251199154849668279D-02];
       

        ZR = z*roots.';
        weights = weights./sqrt(1-roots.*roots);
        zz = exp(1i*z)-1;
        s = sum(weights.*roots);
        sw = sum(weights);
        cvals = (exp(1i*ZR))*weights;
        cvals = sum(cvals,2);
        cvals = cvals - zz*s-sw;
        cvals = 2/pi*(cvals + pi/2+zz);

end

function [cval] = eval_struve_asymp(z)

        zv = -z;
        eps = 1.0d-15;


        cval = 2/pi./zv;
        ct   = cval;

        r1 = abs(ct);

        ifdone = ones(size(zv));

        for i=1:1000
            ct =-ct./zv./zv/4/i/i.* ...
               (2*i*(2*i-1)).^2;
            r2 = abs(ct);

            ifdone = ifdone.*(r2<=r1);
            ifdone = ifdone.*(r2>eps);

            if(sum(ifdone) == 0)
                break
            end
            cval = cval + ct;
            r1 = r2;
        end

        h0 = besselh(0,z);
        indsp = find(real(cval)>0);
        indsn = find(real(cval)<=0);
        cval(indsp) = cval(indsp)-1i*h0(indsp);
        cval(indsn) = cval(indsn)+1i*h0(indsn);

        cval = -1i*cval;
        cval(ifdone==1) = NaN;

end