function cleanup_radios(tx, rx)
    try
        release(tx);
    catch
    end
    try
        release(rx);
    catch
    end
end
